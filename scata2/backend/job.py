from scata2.models import ScataJob


def run_job(pk):
    job = ScataJob.objects.get(pk=pk)
    