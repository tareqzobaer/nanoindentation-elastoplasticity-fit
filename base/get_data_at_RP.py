
from odbAccess import *
import csv

# Open the ODB file
odb_path = 'job_gen_cand_RANDOM.odb'
csv_file = 'job_gen_cand_RANDOM.csv'

odb = openOdb(path=odb_path)

# Access the assembly
assembly = odb.rootAssembly
ref_point_set_name = 'RP_SET'

# Access the set named 'RP_SET'
node_set = assembly.nodeSets[ref_point_set_name]

# Open a CSV file for writing
with open(csv_file, 'wb') as csvfile:
    writer = csv.writer(csvfile)

    # Write the header
    writer.writerow(['Frame', 'Node Label', 'U2', 'RF2'])

    # Assume we are extracting data from the first step
    step = odb.steps['Step-1']

    # Loop through all frames in the step
    for frame in step.frames:
        frame_number = frame.frameId

        # Access the displacement and reaction force field outputs for the current frame
        displacement = frame.fieldOutputs['U']
        reaction_force = frame.fieldOutputs['RF']

        # Loop through the nodes in 'RP_SET'
        for node in node_set.nodes[0]:
            node_label = node.label

            # Get displacement data (U2)
            disp = displacement.getSubset(region=node_set).values[node_label - 1].data
            U2 = disp[1] if len(disp) > 1 else 0.0

            # Get reaction force data (RF2)
            rf = reaction_force.getSubset(region=node_set).values[node_label - 1].data
            RF2 = rf[1] if len(rf) > 1 else 0.0

            # Write data to CSV
            writer.writerow([frame_number, node_label, U2, RF2])

# Close the ODB file
odb.close()
