#!/usr/bin/env python3
with open('Mjav-c.collinearity.kaks.refomatted', 'r') as file:
    lines = file.readlines()


current_group = []
groups_bigger_than = open("groups_bigger_than_10.txt","w")
# Find lines starting with ## Alignment
for line in lines:
    if line.startswith('## Alignment'):
        # With the value of N is bigger than 10, print block
        if len(current_group) > 10:
            for group_line in current_group:
                groups_bigger_than.write(group_line)
            
        current_group = []  # 
    else:
        current_group.append(line)

# 
if len(current_group) > 10:
    for group_line in current_group:
        groups_bigger_than.write(group_line)
        
groups_bigger_than.close()  
