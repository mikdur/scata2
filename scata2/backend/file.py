from scata2.models import ScataFile
import hashlib, sys


def check_file(pk):
    sf = ScataFile.objects.get(pk=pk)
    hash = hashlib.sha256()

    with sf.file.open(mode="rb") as f:
        while True:
            data = f.read(4096)
            if len(data) == 0:
                break
            hash.update(data)

    sf.file_size = sf.file.size / (1024 * 1024)
    sf.sha256 = hash.hexdigest()
    sf.save()
    return True



