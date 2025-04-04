from django.core.files.storage import storages

def get_file_storage():
    if "scata_files" in storages.backends:
        return storages['scata_files']
    else:
        return storages.default_storage
    
def get_work_storage():
    if "scata_work" in storages.backends:
        return storages['scata_work']
    else:
        return storages.default_storage