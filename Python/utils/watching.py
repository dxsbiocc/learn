import pyinotify
import sys


class EventHandler(pyinotify.ProcessEvent):
    
    def my_init(self, file_object=sys.stdout):
        """
        This is your constructor it is automatically called from
        ProcessEvent.__init__(), And extra arguments passed to __init__() would
        be delegated automatically to my_init().
        """
        self._file_object = file_object
        self.sign = False
    
    def process_IN_CREATE(self, event):
        """
        This method processes a specific type of event: IN_CREATE. event
        is an instance of Event.
        """
        if hasattr(event, 'pathname') and event.pathname.endswith('.xls'):
            print('%s have been created.\n' % event.pathname, file=self._file_object)
        
    def process_IN_DELETE(self, event):
        """
        This method processes a specific type of event: IN_DELETE. event
        is an instance of Event.
        """
        print('deleting: %s\n' % event.pathname, file=self._file_object)
        # pass

    def process_IN_CLOSE(self, event):
        """
        This method is called on these events: IN_CLOSE_WRITE and
        IN_CLOSE_NOWRITE.
        """
        # print('closing: %s\n' % event.pathname, file=self._file_object)
        pass
        
    def process_IN_CLOSE_WRITE(self, event):
        """
        This method processes a specific type of event: IN_CLOSE_WRITE
        """
        pass
            
    def process_default(self, event):
        """
        Eventually, this method is called for all others types of events.
        This method can be useful when an action fits all events.
        """
        # print('default processing\n', file=self._file_object)
        pass


def watching(path, exclude_path=None, rec=False, read_freq=0, timeout=None):
    """ watch files or directories
    @args: 
        path:          str or list of str, Path to watch, the path can either be a 
                       file or a directory. Also accepts a sequence (list) of paths.
                       
        exclude_path:  str or list, predicate (boolean function), which returns True 
                       if the current path must be excluded from being watched. This 
                       argument has precedence over exclude_filter passed to the 
                       class' constructor.
                       
        rec:           Recursively add watches from path on all its subdirectories, 
                       set to False by default (doesn't follows symlinks in any case)
                       
        read_freq:     if read_freq == 0, events are read asap, if read_freq is > 0, 
                       this thread sleeps max(0, read_freq - (timeout / 1000)) seconds. 
                       But if timeout is None it may be different because poll is 
                       blocking waiting for something to read.
                       
        timeout:       see read_freq above. If provided, it must be set in milliseconds
    """
    # Instanciate a new WatchManager (will be used to store watches)
    wm = pyinotify.WatchManager()
    # events types
    mask = pyinotify.IN_DELETE | pyinotify.IN_CREATE
    # Associate this WatchManager with a Notifier (will be used to report and process events).
    notifier = pyinotify.Notifier(wm, EventHandler(), read_freq=read_freq, timeout=timeout)
    # Add a new watch on 'path' for some XXX_EVENTS.
    if isinstance(path, str) or isinstance(path, list):
        print("now starting monitor %s." %path)
        wm.add_watch(path, mask, rec=rec, exclude_filter=exclude_path)
    else:
        raise ValueError("the %s seems not valid path" %(path))
    
    # Loop forever and handle events.
    notifier.loop()
            

if __name__ == '__main__':
    path = '~/jupyter/Others/Wechat'
    watching(path)
