class ScataReadsError(Exception):
    """Exception raised for errors when parsing reads
    
    Attributes:
       error - short error code
       message - Human readable error message """

    def __init__(self, error, message="Unknown error"):
        self.error = error
        self.message = message

class ScataFileError(Exception):
    """Exception raised for non-recoverable parse errors
    
    Attributes:
       error - short error code
       message - Human readable error message """
    def __init__(self, error, message="Unknown error"):
        self.error = error
        self.message = message