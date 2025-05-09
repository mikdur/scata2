from django.urls import path, include

from . import views

urlpatterns = [
    path("__reload__/", include("django_browser_reload.urls")),
    path("", views.index, name="scata2-index"),
    path("primers/", views.PrimerListView.as_view(), name="primer-list"),
    path("primers/add/", views.PrimerCreateView.as_view(), name="primer-add"),
    path("primers/<int:pk>/delete/", views.PrimerDeleteView.as_view(), name="primer-delete"),
    path("tagsets/", views.TagSetListView.as_view(), name="tagset-list"),
    path("tagsets/add/", views.TagSetCreateView.as_view(), name="tagset-add"),
    path("tagsets/<int:pk>/delete/", views.TagSetDeleteView.as_view(), name="tagset-delete"),
    path("amplicons/", views.AmpliconListView.as_view(), name="amplicon-list"),
    path("amplicons/add/", views.AmpliconCreateView.as_view(), name="amplicon-add"),
    path("amplicons/<int:pk>/delete/", views.AmpliconDeleteView.as_view(), name="amplicon-delete"),
    path("referenceset/", views.ReferenceSetListView.as_view(), name="referenceset-list"),
    path("referenceset/add/", views.ReferenceSetCreateView.as_view(), name="referenceset-add"),
    path("referenceset/<int:pk>/delete/", views.ReferenceSetDeleteView.as_view(), name="referenceset-delete"),
    path("dataset/", views.DataSetListView.as_view(), name="dataset-list"),
    path("dataset/add/", views.DataSetCreateView.as_view(), name="dataset-add"),
    path("dataset/<int:pk>/delete/", views.DataSetDeleteView.as_view(), name="dataset-delete"),
    path("jobs/", views.JobListView.as_view(), name="job-list"),
    path("jobs/add/", views.JobCreateView.as_view(), name="job-add"),
    path("jobs/<int:pk>/delete/", views.JobDeleteView.as_view(), name="job-delete"),
    path("files/", views.FileListView.as_view(), name="file-list"),
    path("files/add/", views.FileCreateView.as_view(), name="file-add"),
    path("files/<int:pk>/delete/", views.FileDeleteView.as_view(), name="file-delete"),


]

