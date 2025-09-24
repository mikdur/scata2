
import scata2.methods.scata as scata
from scata2.methods.dummy.models import DummyMethod

methods = { "scata" : { "model" : scata.models.ScataMethod,
                        "description": "Scata"},
            "dummy" : { "model" : DummyMethod ,
                        "description" : "Dummy test method" },
}

