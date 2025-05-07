# config.py
import configparser

config = configparser.ConfigParser()
config.read('config.ini')  # assumes it's in the root folder

SECRET_KEY = config.get('REMI', 'SECRET_KEY', fallback='fallback-secret')
NCBI_API_KEY = config.get('REMI', 'NCBI_API_KEY', fallback='')
