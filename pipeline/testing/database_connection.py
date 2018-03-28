def heroku_database():
    import subprocess
    import urllib.parse
    import psycopg2
    result = urllib.parse.urlparse(str(subprocess.check_output("DATABASE_URL=$(heroku config:get DATABASE_URL -a openfoundry) && echo $DATABASE_URL", shell=True), 'utf-8').rstrip())
    connection = psycopg2.connect(
            database = result.path[1:],
            user = result.username,
            password = result.password,
            host =  result.hostname
            )
    return connection

