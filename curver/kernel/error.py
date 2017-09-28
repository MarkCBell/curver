
''' A module for the various different errors that can be raised. '''

class AbortError(Exception):
    ''' An exception for aborting computations with.
    
    This is thrown by clicking 'cancel' on a progress box. '''
    
    def __init__(self, message=None):
        super(AbortError, self).__init__()
        self.message = message
    def __str__(self):
        return str(self.message)

class AssumptionError(Exception):
    ''' An exception for when an assumption is false. '''
    
    def __init__(self, message=None):
        super(AssumptionError, self).__init__()
        self.message = message
    def __str__(self):
        return str(self.message)

