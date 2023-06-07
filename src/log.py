import logging

FORMAT = "[%(asctime)s] [%(levelname)8s] --- %(message)s (%(filename)s:%(lineno)s)"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
LOG_FILE_NAME = 'logs.log'

def logger() -> logging.Logger:
  """
  Returns the current logger.

  :return: Logger.
  :rtype: Logger
  """
  return logging.getLogger()

def initialize_file_log() -> None:
  """
  Initializes the file logger.
  """
  logger = logging.getLogger()
  logger.setLevel(logging.DEBUG)
  if logger.hasHandlers(): 
    logger.handlers = []
  logger.addHandler(_get_handler())

def _get_handler() -> logging.FileHandler:
  """
  Returns the handler for the log file.

  :return: File handler for the log file.
  :rtype: FileHandler
  """
  handler = logging.FileHandler(filename = LOG_FILE_NAME, mode = 'a', encoding = 'utf-8')
  _set_handler_formatter(handler)
  handler.setLevel(logging.INFO)
  return handler 

def _set_handler_formatter(handler: logging.FileHandler) -> None:
  """
  Sets up the formats for the file handler.
  """
  formatter = logging.Formatter(FORMAT, DATE_FORMAT) 
  handler.setFormatter(formatter)

