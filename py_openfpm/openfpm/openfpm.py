import threading
import zmq
import ctypes
import time

requests = {}
context = zmq.Context()
socket = context.socket(zmq.REQ)

def process_completed(message,r_id):
    size_out = int.from_bytes(message[7:7+8], byteorder="little")
    ptr = 7+8
    requests[r_id]["output"] = message[ptr:ptr+size_out].decode()
    ptr += size_out
    size_out = int.from_bytes(message[ptr:ptr + 8], byteorder="little")
    ptr += 8
    requests[r_id]["error"] = message[ptr:ptr + size_out].decode()
    ptr += size_out
    size_out = int.from_bytes(message[ptr:ptr + 8], byteorder="little")
    ptr+= 8
    requests[r_id]["payload"] = message[ptr:ptr+size_out]


def wait_request(r_id):
    completed = False
    while completed is False:
        socket.send(b"queryr:" + ctypes.c_int(r_id))
        message = socket.recv()
        if message[0:7].decode() == "comple:":
            # Process completed message
            process_completed(message,r_id)
            completed = True
        time.sleep(0.1)


def connect(server,port=5555):
    socket.connect("tcp://" + server + ":" + str(port))

def ranks():
    socket.send(b"nranks:")

    message = socket.recv()
    return int.from_bytes(message[7:7+4],byteorder='little')

def grids():
    socket.send(b"lstruc:")
    message = socket.recv()

    ngrids = int.from_bytes(message[7:7+4],byteorder='little')

    ptr = 7+4

    grids = []
    for i in range(0,ngrids):
        grids.append({})
        n_string = int.from_bytes(message[ptr:ptr + 4],byteorder='little')
        ptr += 4
        grids[i]["name"] = message[ptr:ptr + n_string].decode()
        ptr += n_string
        n_string = int.from_bytes(message[ptr:ptr + 4],byteorder='little')
        ptr += 4
        grids[i]["type"] = message[ptr:ptr + n_string].decode()
        ptr += n_string
        grids[i]["dim"] = int.from_bytes(message[ptr:ptr+4],byteorder='little')
        ptr += 4
        n_tot_sizes = int.from_bytes(message[ptr:ptr+4],byteorder='little')
        ptr += 4
        grids[i]["sizes"] = []
        for j in range(0,n_tot_sizes):
            grids[i]["sizes"].append(int.from_bytes(message[ptr:ptr+4],byteorder='little'))
            ptr += 4

    return [d for d in grids if d['name'] == "grid_dist_id"]

def slice(grid, **kwargs):
    

def run_command(code_to_run):

    socket.send(b"runcmd:" + bytes(code_to_run,'utf-8'))

    message = socket.recv()

    if (message[0:7].decode() == "queued:"):
        r_id = int.from_bytes(message[7:7+4],byteorder='little')
        requests[r_id] = {"code": code_to_run, "executed": False}
        wait_request(r_id)
        return requests[r_id]["output"], requests[r_id]["error"]

    return "","Unrecognized answer from the server"


