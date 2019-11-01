"""
helper functions for postgres database connections and transaction
management
"""

import psycopg2
import psycopg2.extras
import os

from CBILDataCommon.Util.FakeSecHead import FakeSecHead #,warning
from ConfigParser import SafeConfigParser

class Database(object):
    '''
    accessor for database connection info
    + database handler
    '''
    def __init__(self, gusConfigFile):
        self.dbh = None # database handler (connection)
        self.dsn = None
        self.user = None
        self.__password = None
        self.__pgpassword = None # placeholder for resetting PGPASSWORD environmental var
        self.dsnConfig = None
        self.connectionString = None

        self.gusConfigFile = os.environ['GUS_HOME'] + "/config/gus.config" \
          if gusConfigFile is None else gusConfigFile

        self.load_database_config()


    def cursor(self, cursorFactory=None):
        '''
        create and return database cursor
        if dictCursor is True, return DictCursor
        '''
        if cursorFactory == 'DictCursor':
            return self.dbh.cursor(cursor_factory=psycopg2.extras.DictCursor)
        if cursorFactory == 'RealDictCursor':
            return self.dbh.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        return self.dbh.cursor()

    def named_cursor(self, name, cursorFactory=None, withhold=True):
        '''
        create and return database cursor
        if dictCursor is True, return DictCursor
        '''
        if cursorFactory == 'DictCursor':
            return self.dbh.cursor(name=name,cursor_factory=psycopg2.extras.DictCursor, withhold=withhold)
        if cursorFactory == 'RealDictCursor':
            return self.dbh.cursor(name=name,cursor_factory=psycopg2.extras.RealDictCursor, withhold=withhold)

        return self.dbh.cursor(name=name, withhold=withhold)


    def set_session(self, readonly=None, autocommit=None):
        '''
        wrapper for setting  one or more parameters for
        the next transactions or statements in the current session.
        pass 'DEFAULT' as the value to reset parameter to server default
        a value of None leaves the setting unchanged
        '''
        self.dbh.set_session(readonly=readonly, autocommit=autocommit)


    def autocommit(self, status=True):
        '''
        sets isolation level (auto-commit mode)
        autocommit must be on for transactions such as create database or vacuum
        default is FALSE
        '''
        self.dbh.autocommit = status


    def set_pgpassword(self):
        '''
        set PGPASSWORD environmental variable 
        b/c postgres does not allow passwords through the command
        line to psql
        '''
        if os.environ.get('PGPASSWORD'):
            self.__pgpassword = os.environ['PGPASSWORD']
        os.environ['PGPASSWORD'] = self.__password

        
    def reset_pgpassword(self):
        '''
        set PGPASSWORD back to original value
        '''
        if self.__pgpassword is not None:
            os.environ['PGPASSWORD'] = self.__pgpassword
        else: # unset pgpassword environmental var
            del os.environ['PGPASSWORD']
            

    def name(self):
        '''
        return database name
        '''
        return self.dsnConfig['dbname']

    
    def port(self):
        '''
        return port
        '''
        if 'port' in self.dsnConfig:
            return self.dsnConfig['port']
        else:
            return None

    def host(self):
        '''
        return database server
        '''
        if 'host' in self.dsnConfig:
            return self.dsnConfig['host']
        else:
            return None

    def load_database_config(self):
        '''
        parse gus config file for DB connection info
        '''
        config_parser = SafeConfigParser()
        config_parser.readfp(FakeSecHead(open(self.gusConfigFile)))

        self.dsn = config_parser.get('section', 'dbiDsn');
        self.user = config_parser.get('section', 'databaseLogin')
        self.__password = config_parser.get('section', 'databasePassword')

        self.dsn = self.dsn.replace('dbi:Pg:', '')
        self.dsn = self.dsn.replace('DBI:Pg:', '')
        self.dsnConfig = dict(param.split("=") for param in self.dsn.split(";"))


    def connect(self):
        '''
        create database connection
        '''
        self.connectionString = "user='" + self.user + "'";
        self.connectionString = self.connectionString + " password='" + self.__password + "'"
        self.connectionString = self.connectionString + " dbname='" + self.dsnConfig['dbname'] + "'"
        if 'host' in self.dsnConfig:
            self.connectionString = self.connectionString + " host='" + self.dsnConfig['host'] + "'"
        if 'port' in self.dsnConfig:
            self.connectionString = self.connectionString + " port='" + self.dsnConfig['port'] + "'"
    
        self.dbh = psycopg2.connect(self.connectionString)


    def close(self):
        '''
        close the database connection
        '''
        self.dbh.close()


    def rollback(self):
        '''
        rollback any changes
        '''
        self.dbh.rollback()


    def commit(self):
        '''
        commit any changes
        '''
        self.dbh.commit()


  
