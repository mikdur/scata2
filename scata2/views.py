from django.shortcuts import render
from django.urls import reverse_lazy
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.views.generic import ListView
from django.views.generic.edit import CreateView, DeleteView
from scata2.models import ScataPrimer


class OwnerCheckMixin(LoginRequiredMixin, UserPassesTestMixin):
    def test_func(self):
        return self.get_object().owner == self.request.user


# Root page
@login_required
def index(request):
    return render(request, "scata2/index.html", {})

#################################
#  Primer views
#################################

class PrimerListView(LoginRequiredMixin,ListView):
    model = ScataPrimer

    def get_queryset(self):
        return ScataPrimer.objects.filter(owner=self.request.user)

class PrimerCreateView(LoginRequiredMixin,CreateView):
    model = ScataPrimer
    fields = ["short_name", "sequence", "mismatches", "description"]

    def form_valid(self, form):
        form.instance.owner = self.request.user
        return super().form_valid(form)

class PrimerDeleteView(OwnerCheckMixin,DeleteView):
    model = ScataPrimer
    success_url = reverse_lazy("primer-list")
