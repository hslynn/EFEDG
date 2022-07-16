import matplotlib.pyplot as plt

def read_data(filename, value_idx, max_time):
    time_list = []
    value_list = []
    
    with open(filename,"r") as f:
        for line in f.readlines():
            obj_list = line.split(" ")
            if float(obj_list[0]) <= max_time:
                time_list.append(float(obj_list[0]))
                value_list.append(float(obj_list[value_idx]))
    return [time_list, value_list]
