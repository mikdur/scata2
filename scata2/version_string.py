import Bio, numpy, django, scata2, sklearn

def package_versions():
    return "scata2={sc}, biopython={bp}, numpy={np}, scikit-learn={sk}, django={dj}".format(
        bp=Bio.__version__, 
        np=numpy.__version__,
        sk=sklearn.__version__,
        dj=django.__version__,
        sc=scata2.__version__
    )