# Copy this file to config.py and make changes

BASE_PATH = "/Users/conarymeyer/bionet-synthesis"

# SQL_USERNAME = "openfoundry-dev"
# SQL_PASSWORD = "founDRY"
SQL_USERNAME = "openfoundry"
SQL_PASSWORD = "freegenestomakegenesfree"

# CONNECTION_STRING = 'postgresql+psycopg2://{}:{}@freegenes-openfoundry.cwtlxuukykrr.us-east-1.rds.amazonaws.com:5432/openfoundry'.format(SQL_USERNAME,SQL_PASSWORD)
CONNECTION_STRING = 'sqlite:///template.db'
# CONNECTION_STRING = 'sqlite:///:memory:'

BIONET_RPC_URL = "http://localhost:8000/rpc"
BIONET_USERNAME = None
BIONET_PASSWORD = None

LOGGING_CONFIG = {
    'version': 1,
    'formatters': {
      'simple': {
        'format': '%(asctime)s %(filename)s %(name)s %(levelname)s %(message)s'
      }
    },
    'handlers': {
      'console': {
        'class': 'logging.StreamHandler',
        'level': 'DEBUG',
        'formatter': 'simple',
        'stream': 'ext://sys.stdout',
      },
      'file': {
        'class': 'logging.FileHandler',
        'level': 'DEBUG',
        'formatter': 'simple',
        'filename': BASE_PATH + '/log/pipeline.log'
      }
    },
    'loggers': {
      '__main__': {
        'level': 'DEBUG',
        'handlers': ['console', 'file'],
        'propagate': False
      }
    },
    'root': {
      'level': 'DEBUG',
      'handlers': ['console']
    }
}
