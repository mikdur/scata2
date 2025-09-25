
import scata2.methods.scata as scata
import scata2.methods.dummy as dummy

methods = {"scata": {"model": scata.models.ScataMethod,
                     "form": scata.models.ScataMethodForm,
                     "description": "Scata"},
           "dummy": {"model": dummy.models.DummyMethod,
                     "form": dummy.models.DummyMethodForm,
                     "description": "Dummy test method"},
           }
