
import scata2.methods.scata.models as scata_models
import scata2.methods.dummy.models as dummy_models

methods = {"scata": {"model": scata_models.ScataMethod,
                     "form": scata_models.ScataMethodForm,
                     "description": "Scata classic"},
           "dummy": {"model": dummy_models.DummyMethod,
                     "form": dummy_models.DummyMethodForm,
                     "description": "Dummy test"},
           }
