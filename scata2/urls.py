from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name="scata2-index"),
    path("primers/", views.PrimerListView.as_view(), name="primer-list"),
    path("primers/add/", views.PrimerCreateView.as_view(), name="primer-add"),
    path("primers/<int:pk>/delete/", views.PrimerDeleteView.as_view(), name="primer-delete"),
    path("tagset/", views.TagSetListView.as_view(), name="tagset-list"),
    path("tagset/add/", views.TagSetCreateView.as_view(), name="tagset-add"),
    path("tagset/<int:pk>/delete/", views.TagSetDeleteView.as_view(), name="tagset-delete"),
    path("files/", views.FileListView.as_view(), name="file-list"),
    path("files/add/", views.FileCreateView.as_view(), name="file-add"),
    path("files/<int:pk>/delete/", views.FileDeleteView.as_view(), name="file-delete"),


]

