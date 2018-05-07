
from jsonrpc_requests import Server, ProtocolError
import json

class Bionet():

    # update the status in the bionet

    def __init__(self, url, debug=None):

        self.debug = debug
        
        # bionet nodes have a json_rpc endpoint at /rpc
        self.bionet = Server(url)    

    def login(self, username, password):

        # the authentication token is automatically stored in a cookie
        try:
            res = self.bionet.login(username, password)
            if self.debug:
                print("[bionet] login successful")
        except:
            print("[bionet] login failed", file=sys.stderr)
            return False

        return True

    def update_status(self, virtual_id, status_str):

        try:
            res = self.bionet.get(virtual_id)
            if self.debug:
                print("[bionet] found biomaterial with id: %s" % virtual_id)
        
            material = json.loads(res[0])
        except:
            print("[bionet] failed to get biomaterial with id: %s" % virtual_id, file=sys.stderr)
            return False
        
        try:
            material['freegenes']['status']['current_status'] = status_str
            
            self.bionet.saveVirtual(material)
        
        except:
            print("[bionet] failed to get biomaterial with id: %s" % virtual_id, file=sys.stderr)
            return False

        return True




        
