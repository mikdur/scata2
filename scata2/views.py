from django.shortcuts import render
from django.urls import reverse_lazy
from django.http import HttpResponseRedirect
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.views.generic import ListView
from django.views.generic.edit import CreateView, DeleteView
from scata2.models import ScataPrimer, ScataTagSet

import sys


# Mixin to ensure user is logged in and that the user is the owner of
# the object.

class OwnerCheckMixin(UserPassesTestMixin):
    def test_func(self):
        o = self.get_object()
        return o.owner == self.request.user and not o.deleted

class ListOwnedView(ListView,LoginRequiredMixin):
    def get_queryset(self):
        return self.model.objects.filter(owner=self.request.user,
                                         deleted=False)

class DeleteToTrashView(DeleteView,LoginRequiredMixin,
                       OwnerCheckMixin):
    def form_valid(self, form):
        self.object.deleted=True
        self.object.save()
        success_url = self.get_success_url()
        return HttpResponseRedirect(success_url)
    

# Root page
@login_required
def index(request):
    return render(request, "scata2/index.html", {})

#################################
#  Primer views
#################################

class PrimerListView(ListOwnedView):
    model = ScataPrimer

    
class PrimerCreateView(CreateView,LoginRequiredMixin):
    model = ScataPrimer
    fields = ["short_name", "sequence", "mismatches", "description"]

    def form_valid(self, form):
        form.instance.owner = self.request.user
        return super().form_valid(form)

class PrimerDeleteView(DeleteToTrashView):
    model = ScataPrimer
    success_url = reverse_lazy("primer-list")

    


#################################
#  Tagset views
#################################

class TagSetListView(ListOwnedView):
    model = ScataTagSet

class TagSetCreateView(CreateView,LoginRequiredMixin):
    model = ScataTagSet
    fields = ["name", "tagset_file"]

    def form_valid(self, form):
        form.instance.owner = self.request.user
        return super().form_valid(form)
    

class TagSetDeleteView(DeleteToTrashView):
    model = ScataTagSet
    success_url = reverse_lazy("tagset-list")

    
