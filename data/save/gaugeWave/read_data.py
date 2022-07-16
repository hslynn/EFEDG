import matplotlib.pyplot as plt

def read_data(filename, value_idx, max_time):
    time_list = []
    value_list = []
    
    with open(filename,"r") as f:
        for line in f.readlines():
            obj_list = line.split(" ")
            time = float(obj_list[0]);
            if time <= max_time:
                if len(time_list) == 0 or time_list[len(time_list)-1] <= time:
                    time_list.append(time)
                    value_list.append(float(obj_list[value_idx]))
    return [time_list, value_list]
