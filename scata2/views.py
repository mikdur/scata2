from django.shortcuts import render
from django.urls import reverse_lazy
from django.http import HttpResponseRedirect, HttpResponse, JsonResponse
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.views.generic import ListView, DetailView
from django.views.generic.edit import CreateView, DeleteView
from django.db.models import Q
from django.core.exceptions import SuspiciousOperation
from scata2.models import ScataFile, ScataPrimer, ScataTagSet, ScataAmplicon, \
                          ScataReferenceSet, ScataRefsetErrorType, \
                          ScataDataset, ScataErrorType, ScataTagStat, \
                          ScataJob, ScataModel
import scata2.backend
from scata2.methods import methods as clustering_methods

import django_q.tasks as q2
import csv, urllib


class FilteredCreateView(LoginRequiredMixin, CreateView):
    # Set owner on save
    def form_valid(self, form):
        form.instance.owner = self.request.user
        return super().form_valid(form)

    # Filter any field querysets subclassed from ScataModel to only include owned
    # objects
    def get_form(self, *args, **kwargs):
        form = super().get_form(*args, **kwargs)  # Get the form as usual
        user = self.request.user
        for k in form.fields.keys():
            if ( hasattr(form.fields[k], "queryset") and 
                issubclass(form.fields[k].queryset.model, ScataModel) ):
                form.fields[k].queryset = \
                    form.fields[k].queryset.filter(Q(owner=user, deleted=False) |
                                                   Q(public=True, deleted=False))
        return form


class ListOwnedView(LoginRequiredMixin, ListView):
    def get_queryset(self):
        return self.model.objects.filter(Q(owner=self.request.user,
                                           deleted=False) |
                                         Q(public=True, deleted=False))

class DeleteToTrashView(LoginRequiredMixin,
                       UserPassesTestMixin,DeleteView):
    def test_func(self):
        o = self.get_object()
        return o.owner == self.request.user and not o.deleted
    
    def form_valid(self, form):
        self.object.deleted=True
        self.object.save()
        success_url = self.get_success_url()
        return HttpResponseRedirect(success_url)

class OwnedDetailView(LoginRequiredMixin, UserPassesTestMixin, 
                      DetailView):
    def test_func(self):
        o = self.get_object()
        return o.owner == self.request.user and not o.deleted

# Mixin to render JSON response
class JSONResponseMixin:

    def render_to_response(self, context, **response_kwargs):
        return JsonResponse(self.get_data(context), **response_kwargs)

    # To be overidden
    def get_data(self, context):
        return context

# Mixin to render CSV response
# Overide get_data with a function that returns
# a list of lists.
class CSVResponseMixin:

    def render_to_response(self, context, **response_kwargs):
        response = HttpResponse(content_type="text/csv",
            headers={"Content-Disposition": 'attachment; filename="{fn}"'.format(fn=self.get_filename(context))},
            **response_kwargs)
        writer=csv.writer(response)
        for row in self.get_data(context):
            writer.writerow(row)
        return response

    # To be overidden
    def get_data(self, context):
        return ["",""]
    
    def get_filename(self, context):
        return "data.csv";


# Root page
@login_required
def index(request):
    return render(request, "scata2/index.html", {})


#################################
#  File views
#################################

class FileListView(ListOwnedView):
    model = ScataFile

class FileCreateView(FilteredCreateView):
    model = ScataFile
    fields = ["name", "file", "description"]

    def form_valid(self, form):
        # Call the parent's form_valid() to save the form
        response = super().form_valid(form)
        task_id = q2.async_task(scata2.backend.check_file, self.object.pk,
                      task_name="file checksum pk={id}".format(id=self.object.pk))
        return response

class FileDeleteView(DeleteToTrashView):
    model = ScataFile
    success_url = reverse_lazy("file-list")


#################################
#  Primer views
#################################

class PrimerListView(ListOwnedView):
    model = ScataPrimer

    
class PrimerCreateView(FilteredCreateView):
    model = ScataPrimer
    fields = ["name", "description", "sequence", "mismatches", "description"]

class PrimerDeleteView(DeleteToTrashView):
    model = ScataPrimer
    success_url = reverse_lazy("primer-list")

    


#################################
#  Tagset views
#################################

class TagSetListView(ListOwnedView):
    model = ScataTagSet

class TagSetCreateView(FilteredCreateView):
    model = ScataTagSet
    fields = ["name", "description", "tagset_file"]

    def form_valid(self, form):
        # Call the parent's form_valid() to save the form
        response = super().form_valid(form)
        task_id = q2.async_task(scata2.backend.parse_tagset, self.object.pk,
                      task_name="tagset pk={id}".format(id=self.object.pk))
        return response

    
class TagSetDeleteView(DeleteToTrashView):
    model = ScataTagSet
    success_url = reverse_lazy("tagset-list")


#################################
#  Amplicon views
#################################

class AmpliconListView(ListOwnedView):
    model = ScataAmplicon

class AmpliconCreateView(FilteredCreateView):
    model = ScataAmplicon
    fields = ["name", "description", "five_prime_primer", "five_prime_tag",
              "three_prime_primer", "three_prime_tag", 
              "min_length", "max_length"]
    
    def get_form(self, *args, **kwargs):
        form = super().get_form(*args, **kwargs)  # Get the form as usual
        for k in ['five_prime_tag', 'three_prime_tag']:
            form.fields[k].queryset = \
                form.fields[k].queryset.filter(validated=True)
        return form

class AmpliconDeleteView(DeleteToTrashView):
    model = ScataAmplicon
    success_url = reverse_lazy("tagset-list")

#################################
#  RefSeq views
#################################

class ReferenceSetListView(ListOwnedView):
    model = ScataReferenceSet

class ReferenceSetCreateView(FilteredCreateView):
    model = ScataReferenceSet
    fields = ["name", "description", "refseq_file", "amplicon"]

    def form_valid(self, form):
        # Call the parent's form_valid() to save the form
        response = super().form_valid(form)
        task_id = q2.async_task(scata2.backend.check_refset, self.object.pk,
                      task_name="refset pk={id}".format(id=self.object.pk))
        return response
    
class ReferenceSetDeleteView(DeleteToTrashView):
    model = ScataReferenceSet
    success_url = reverse_lazy("referenceset-list")


class ReferenceSetDetailView(OwnedDetailView):
    model = ScataReferenceSet

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['errors'] = ScataRefsetErrorType.objects.all().filter(refset=context['object'].pk)
        print(context['errors'])
        return context




#################################
#  Dataset views
#################################

class DataSetListView(ListOwnedView):
    model = ScataDataset

class DataSetCreateView(FilteredCreateView):
    model = ScataDataset
    fields = ["name", "description", "amplicon", "min_qual", "mean_qual",
              "filter_method", "file_types", "file1", "file2",
              "kmer_size", "kmer_hsp_count", "kmer_shared"]
    
    def form_valid(self, form):
        # Call the parent's form_valid() to save the form
        response = super().form_valid(form)
        task_id = q2.async_task(scata2.backend.check_dataset, self.object.pk,
                      task_name="dataset pk={id}".format(id=self.object.pk))
        return response
    
class DataSetDeleteView(DeleteToTrashView):
    model = ScataDataset
    success_url = reverse_lazy("dataset-list")

class DataSetDetailView(OwnedDetailView):
    model = ScataDataset

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['errors'] = ScataErrorType.objects.all().filter(dataset=context['object'].pk)
        context['tags'] = ScataTagStat.objects.all().filter(dataset=context['object'].pk)
        context['discarded'] = context['object'].seq_total - context['object'].seq_count
        return context

class DataSetTagsCSVView(CSVResponseMixin, DataSetDetailView):

    def get_data(self, context):
        data = context['tags'].order_by("count").values().values()
        columns = s = ['tag', 'count'] + \
            sorted(list(set(data.first().keys()) - 
                        set(['id', 'dataset_id','tag', 'count',
                             "min_gc", "max_gc"])))
        return [columns] + \
            [[q[b] for b in columns] for q  in data ]
    
    def get_filename(self, context):
        o = context['object']
        return urllib.parse.quote("dataset_{i}_{n}.csv".format(i=o.id, n=o.name),
                                  safe="")
    
class DataSetTagsJSONView(JSONResponseMixin, DataSetDetailView):

    def get_data(self, context):
        return dict(data=list(context['tags'].order_by("count").values()))

class DataSetTagsPCAJSONView(JSONResponseMixin, DataSetDetailView):

    def get_data(self, context):
        data = list(context['tags'].filter(in_pca=True).order_by("count").values())
        data = [a | {"rev_freq":a['reversed'] / a['count']} for a in data]
        return dict(data=data)

#################################
#  Job views
#################################

class JobListView(ListOwnedView):
    model = ScataJob

class JobCreateView(FilteredCreateView):
    model = ScataJob
    fields = ["name", "description", "datasets", 
              "refsets", "repseqs", "method",]

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        context['method_forms'] = {k: {'form':v["form"](),
                                       'description':v['description']} 
                                       for k, v in clustering_methods.items()}
        return context

    # Catch class based form validation to take care of method-specific 
    # form data and validate.

    def form_valid(self, form):
        form.instance.owner = self.request.user
        method = form.data['method']
        if method not in clustering_methods:
            raise SuspiciousOperation("Unknown scata method")


        method_form = clustering_methods[method]['form'](self.request.POST)
        if method_form.errors:
            return super().form_invalid(form)

        self.object = form.save()
        method_form.instance.job = self.object
        method_form.save()

        q2.async_task(scata2.backend.run_job, self.object.pk,
                      task_name="job pk={id}".format(id=self.object.pk))

        return HttpResponseRedirect(self.get_success_url())


class JobDetailView(OwnedDetailView):
    model = ScataJob

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        m=self.object.method
        context['method_object']=clustering_methods[m]['model'].objects.all().\
            filter(job=self.object.pk)[0]
        return context

class JobDeleteView(DeleteToTrashView):
    model = ScataJob
    success_url = reverse_lazy("job-list")
