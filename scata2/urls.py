from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name="index"),
    path("primers/", views.PrimerListView.as_view(), name="primer-list"),
    path("primers/add/", views.PrimerCreateView.as_view(), name="primer-add"),
    path("primers/<int:pk>/delete/", views.PrimerDeleteView.as_view(), name="primer-delete")
]

