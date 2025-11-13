
import scata2.methods.scata.models as scata_models
import scata2.methods.dummy.models as dummy_models

methods = {"scata": {"model": scata_models.ScataScataMethod,
                     "form": scata_models.ScataScataMethodForm,
                     "description": "Scata classic"},
           "dummy": {"model": dummy_models.ScataDummyMethod,
                     "form": dummy_models.ScataDummyMethodForm,
                     "description": "Dummy test"},
           }
