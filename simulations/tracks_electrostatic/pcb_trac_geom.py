import gmsh, sys

MAX_ELEMENT_SIZE    = .2
GROUND_PLANE_WIDTH  = 10
TRACKS_WIDTH        = 3
AIR_BOX_HEIGHT      = 6
EXTRUDE_LENGTH      = 3
##########################################################################
gmsh.initialize()
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.set_number("Mesh.CharacteristicLengthMax", MAX_ELEMENT_SIZE)
model = gmsh.model
model.add("PCB")
kernel = gmsh.model.occ
##########################################################################



lower_track = kernel.add_rectangle(-TRACKS_WIDTH/2, 0, 0, TRACKS_WIDTH, 1)
upper_track = kernel.add_rectangle(-TRACKS_WIDTH/2, 1.5, 0, TRACKS_WIDTH, 1)
ground_plane = kernel.add_rectangle(-GROUND_PLANE_WIDTH/2, -1.5, 0, GROUND_PLANE_WIDTH, 1)
air_box = kernel.add_rectangle(-GROUND_PLANE_WIDTH/2, -1.5, 0, GROUND_PLANE_WIDTH, AIR_BOX_HEIGHT)

print(f"Fragmenting {(lower_track, upper_track, ground_plane, air_box)}")
outDimTags, outDimTagsMap = kernel.fragment(
    list(zip(
        (2,2,2,2),
        (lower_track, upper_track, ground_plane, air_box)
    )),
    []
)

print(f"{outDimTags = }")
print(f"{outDimTagsMap = }")

extruded_dim_tags = kernel.extrude(outDimTags, 0, 0, EXTRUDE_LENGTH)

kernel.synchronize()

model.add_physical_group(3, [lower_track], tag=1,  name="lower_track")
model.add_physical_group(3, [upper_track], tag=2, name="upper_track")
model.add_physical_group(3, [ground_plane], tag=3, name="ground_plane")
model.add_physical_group(3, [air_box], tag=4, name="air_box")


model.add_physical_group(2, [14], tag=5, name="anode")
model.add_physical_group(2, [2], tag=6, name="cathode")
model.add_physical_group(2, [15], tag=7, name="gnd")

all_volumes = [lower_track, upper_track, ground_plane, air_box]
model.add_physical_group(3, all_volumes, tag=8, name="all")

##########################################################################
model.mesh.field.set_number
model.mesh.generate()
gmsh.write("simulations/tracks_electrostatic/pcb_track.msh")
