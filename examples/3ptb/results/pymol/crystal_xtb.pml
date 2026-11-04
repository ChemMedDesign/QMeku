# QMeku -> PyMOL | pose crystal, xtb_pair (kcal/mol), scale -13.14 .. 0 .. +1.51
# red = most negative, light grey ~ 0, light blue = positive
# Run from this folder:  pymol crystal_xtb.pml
reinitialize
load receptor.pdb, receptor
load ligand_crystal.sdf, ligand
hide everything
bg_color white
show cartoon, receptor and polymer
color grey90, receptor
set cartoon_transparency, 0.7, receptor
set cartoon_fancy_helices, 1
show sticks, ligand
show spheres, ligand
color paleyellow, ligand
util.cnc ligand
set sphere_scale, 0.2, ligand
set sphere_scale, 0.1, ligand and hydro
set stick_radius, 0.14, receptor
set stick_radius, 0.18, ligand
set cartoon_loop_radius, 0.12
set cartoon_oval_length, 0.9
set cartoon_rect_length, 0.9
set_color qm_000, [0.900, 0.000, 0.000]
show sticks, (receptor and polymer and resi 189)
util.cnc (receptor and polymer and resi 189)
color qm_000, (receptor and polymer and resi 189) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 189)
set_color qm_001, [0.871, 0.498, 0.498]
show sticks, (receptor and polymer and resi 190)
util.cnc (receptor and polymer and resi 190)
color qm_001, (receptor and polymer and resi 190) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 190)
set_color qm_002, [0.865, 0.597, 0.597]
show sticks, (receptor and polymer and resi 191)
util.cnc (receptor and polymer and resi 191)
color qm_002, (receptor and polymer and resi 191) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 191)
set_color qm_003, [0.853, 0.804, 0.804]
show sticks, (receptor and polymer and resi 192)
util.cnc (receptor and polymer and resi 192)
color qm_003, (receptor and polymer and resi 192) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 192)
set_color qm_004, [0.862, 0.642, 0.642]
show sticks, (receptor and polymer and resi 195)
util.cnc (receptor and polymer and resi 195)
color qm_004, (receptor and polymer and resi 195) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 195)
set_color qm_005, [0.865, 0.588, 0.588]
show sticks, (receptor and polymer and resi 213)
util.cnc (receptor and polymer and resi 213)
color qm_005, (receptor and polymer and resi 213) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 213)
set_color qm_006, [0.855, 0.773, 0.773]
show sticks, (receptor and polymer and resi 214)
util.cnc (receptor and polymer and resi 214)
color qm_006, (receptor and polymer and resi 214) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 214)
set_color qm_007, [0.868, 0.545, 0.545]
show sticks, (receptor and polymer and resi 215)
util.cnc (receptor and polymer and resi 215)
color qm_007, (receptor and polymer and resi 215) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 215)
set_color qm_008, [0.866, 0.573, 0.573]
show sticks, (receptor and polymer and resi 216)
util.cnc (receptor and polymer and resi 216)
color qm_008, (receptor and polymer and resi 216) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 216)
set_color qm_009, [0.874, 0.435, 0.435]
show sticks, (receptor and polymer and resi 219)
util.cnc (receptor and polymer and resi 219)
color qm_009, (receptor and polymer and resi 219) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 219)
set_color qm_010, [0.877, 0.385, 0.385]
show sticks, (receptor and polymer and resi 220)
util.cnc (receptor and polymer and resi 220)
color qm_010, (receptor and polymer and resi 220) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 220)
set_color qm_011, [0.567, 0.784, 0.991]
show sticks, (receptor and polymer and resi 221A)
util.cnc (receptor and polymer and resi 221A)
color qm_011, (receptor and polymer and resi 221A) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 221A)
set_color qm_012, [0.764, 0.830, 0.893]
show sticks, (receptor and polymer and resi 226)
util.cnc (receptor and polymer and resi 226)
color qm_012, (receptor and polymer and resi 226) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 226)
set_color qm_013, [0.620, 0.796, 0.965]
show sticks, (receptor and polymer and resi 227)
util.cnc (receptor and polymer and resi 227)
color qm_013, (receptor and polymer and resi 227) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 227)
set_color qm_014, [0.550, 0.780, 1.000]
show sticks, (receptor and polymer and resi 228)
util.cnc (receptor and polymer and resi 228)
color qm_014, (receptor and polymer and resi 228) and elem C
set cartoon_transparency, 0.0, (receptor and polymer and resi 228)
alter receptor, b=0.0
alter (receptor and polymer and resi 189), b=-13.1383
alter (receptor and polymer and resi 190), b=-3.0215
alter (receptor and polymer and resi 191), b=-1.7412
alter (receptor and polymer and resi 192), b=-0.1003
alter (receptor and polymer and resi 195), b=-1.2596
alter (receptor and polymer and resi 213), b=-1.8512
alter (receptor and polymer and resi 214), b=-0.2401
alter (receptor and polymer and resi 215), b=-2.3750
alter (receptor and polymer and resi 216), b=-2.0256
alter (receptor and polymer and resi 219), b=-3.9840
alter (receptor and polymer and resi 220), b=-4.8011
alter (receptor and polymer and resi 221A), b=1.3692
alter (receptor and polymer and resi 226), b=0.1876
alter (receptor and polymer and resi 227), b=0.9678
alter (receptor and polymer and resi 228), b=1.5108
alter ligand, b=0.0
set label_font_id, 7
set label_size, 13
set label_color, black
set label_outline_color, white
set label_position, [1.5, 1.5, 4]
label ((receptor and polymer and resi 189) and name CA), "Asp189"
distance qm_polar_1, ligand and index 9, (receptor and polymer and resi 189) and name OD2
label ((receptor and polymer and resi 220) and name CA), "Cyx220"
label ((receptor and polymer and resi 219) and name CA), "Gly219"
distance qm_polar_2, ligand and index 8, (receptor and polymer and resi 219) and name O
label ((receptor and polymer and resi 190) and name CA), "Ser190"
distance qm_polar_3, ligand and index 9, (receptor and polymer and resi 190) and name OG
label ((receptor and polymer and resi 215) and name CA), "Trp215"
label ((receptor and polymer and resi 216) and name CA), "Gly216"
label ((receptor and polymer and resi 213) and name CA), "Val213"
label ((receptor and polymer and resi 191) and name CA), "Cyx191"
set dash_gap, 0.25
set dash_radius, 0.07
set dash_color, black
set dash_round_ends, 0
hide labels, qm_polar_*
select qm_pocket, ligand or (receptor and polymer and resi 189) or (receptor and polymer and resi 190) or (receptor and polymer and resi 191) or (receptor and polymer and resi 192) or (receptor and polymer and resi 195) or (receptor and polymer and resi 213) or (receptor and polymer and resi 214) or (receptor and polymer and resi 215) or (receptor and polymer and resi 216) or (receptor and polymer and resi 219) or (receptor and polymer and resi 220) or (receptor and polymer and resi 221A) or (receptor and polymer and resi 226) or (receptor and polymer and resi 227) or (receptor and polymer and resi 228)
orient qm_pocket
zoom qm_pocket, 4
set ray_shadows, 0
set antialias, 2
set depth_cue, 0
set specular, 0.2
hide everything, (receptor or ligand) and hydro and not neighbor elem N+O+S
deselect
