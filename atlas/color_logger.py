import os, sys
import logging, traceback

# root logger
logger = logging.getLogger()

# logger= logging


grey = "\x1b[38;21m"
green = "\x1b[32;21m"
yellow = "\x1b[33;21m"
red = "\x1b[31;21m"
bold_red = "\x1b[31;1m"
reset = "\x1b[0m"

prefix = "[Atlas] "


class ColorFormatter(logging.Formatter):
    def __init__(
        self,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)",
    ):

        self.FORMATS = {
            logging.DEBUG: prefix + grey + format + reset,
            logging.INFO: prefix + green + format + reset,
            logging.WARNING: prefix + yellow + format + reset,
            logging.ERROR: prefix + red + format + reset,
            logging.CRITICAL: prefix + red + format + reset,
        }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


#
logging_format = "%(levelname)s: %(message)s"
# datefmt="%Y-%m-%d %H:%M"
#
#
#
#
#
# fileHandler = logging.FileHandler("atlas.log",mode='w')
# fileHandler.setFormatter(logging.Formatter(logging_format))
# fileHandler.setLevel(logging.DEBUG)
#
#
#
# creat console logging
consoleHandler = logging.StreamHandler()
consoleHandler.setLevel(logging.INFO)
consoleHandler.setFormatter(ColorFormatter(logging_format))


#

## Define logging
logging.basicConfig(
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M",
    format=logging_format,
    handlers=[consoleHandler],
)
logging.captureWarnings(True)


# create logging for atlas


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception
