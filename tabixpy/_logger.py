import logging

from ._logger import *

# Custom formatter
#https://stackoverflow.com/questions/14844970/modifying-logging-message-format-based-on-message-logging-level-in-python3
class CustomConsoleFormatter(logging.Formatter):
    _fmts = {
        # logging.NOSET   : "NOSET   : %(msg)s",
        logging.DEBUG   : "DEBUG   : %(asctime)8s - %(name)-7s/%(module)-7s/%(funcName)-12s/%(lineno)3d - %(msg)s",
        logging.INFO    : "INFO    : %(msg)s",
        logging.WARNING : "WARNING : %(msg)s",
        logging.ERROR   : "ERROR   : %(name)-7s/%(module)-7ss/%(funcName)-12s/%(lineno)3d - %(msg)s",
        logging.CRITICAL: "CRITICAL: %(name)-7s/%(module)-7ss/%(funcName)-12s/%(lineno)3d - %(msg)s",
    }
    # format  = "%(asctime)8s - %(name)-7s/%(funcName)-12s/%(lineno)3d - %(levelname)-7s - %(message)s",
    # datefmt = "%H:%M:%S",
    # datefmt='%Y-%m-%d %H:%M:%S',

    def __init__(self, fmt="%(levelno)d: %(msg)s", datefmt="%H:%M:%S"):
        super().__init__(fmt=fmt, datefmt=datefmt, style='%')

    def format(self, record):
        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        self._style._fmt = self._fmts.get(record.levelno, format_orig)

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result

def setLogLevel(level):
    logger.setLevel(level)

def getLogLevel():
    return logging.getLevelName(logger.getEffectiveLevel())



# Set up a logger
logger          = logging.getLogger('tabixpy')

my_formatter    = CustomConsoleFormatter(datefmt="%H:%M:%S",)

console_handler = logging.StreamHandler()
console_handler.setFormatter(my_formatter)

logger.addHandler(console_handler)
logger.setLevel(logging.INFO)
