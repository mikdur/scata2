from django.db import models


class Method(models.Model):
    job = models.OneToOneField("scata2.ScataJob", on_delete=models.CASCADE,
                               related_name="job_method")

    @classmethod
    def run_job(cls, job):
        instance = cls.objects.get(job=job)
        instance.cluster()
