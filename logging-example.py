import logging
import logging.config
from config import *

print(BASE)
print(LOGGING_CONFIG)

logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.info("build-01 BBF10K_00001 message"')
