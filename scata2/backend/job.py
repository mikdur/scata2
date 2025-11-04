from scata2.models import ScataJob
from scata2.methods import methods as clustering_methods


def run_job(pk):
    job = ScataJob.objects.get(pk=pk)

    model = clustering_methods[job.method]['model']

    model.run_job(job=job.pk)
