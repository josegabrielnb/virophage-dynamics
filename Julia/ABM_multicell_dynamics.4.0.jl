### Agent-based model of cell/virus/virophage system and multicellularity

### Parameters set for a population of single-celled organisms and neutral virophages.
### To simulate other scenarios, change the f, cell_number and programmed_cell_death 
### variables accordingly.

### You may need to install the necessary Julia packages. The option to produce an
### animated gif is turned off, you may modify the code to activate this functionality.

### Developed by Jose Gabriel Nino Barreat, jose.ninobarreat@zoo.ox.ac.uk (2022)

## Packages

using Plots
using Random
using Distributions: Uniform
using Distributions: Binomial
using NearestNeighbors
using DelimitedFiles
using Distances
using StatsBase

#Random.seed!(1234)

## Agents

# Create agent types

#define cell agents
mutable struct cell
    t::Int64 #time of birth
    o::Int64 #organism identifier
    s::Array{Int64,1} #state of cell from time t'
    x::Float64
    y::Float64
    z::Float64
end

#=
Cell states:
1 --> uninfected cell
2 --> infected by giant virus
3 --> infected by exogenous virophage
4 --> virophage genome integration
5 --> infected by giant virus and exogenous virophage (inhibition)
6 --> infected by giant virus and endogenous virophage
7 --> infected by endogenous and exogenous virophage
8 --> infected by endogenous, exogenous virophages and giant virus (inhibition)
=#

#define virus agents
mutable struct virus
    t::Int64 #time of birth
    x::Float64
    y::Float64
    z::Float64
end

#define virophage agents
mutable struct virophage
    t::Int64 #time of birth
    x::Float64
    y::Float64
    z::Float64
end

## Parameters

# Model parameters

N0 = 1024 #initial population size
cell_radius = 20 #cell radius in micrometers
virus_radius = 0.25 # virus radius in micrometers
virophage_radius = 0.04 # virophage radius in micrometers
space = 1e4 #space in one dimension in micrometers
V0 = 512 #initial virus population size
M0 = 2048 #initial virophage population size
k1 = 1e-7 #rate of virus decay
k2 = 8e-8 #rate of virophage decay
integration_p = 0.1 #0.1 #integration probability
f = 1 #virophage inhibition of giant virus replication (1:no inhibition, 0.01:inhibition)

#travelled distances of cells and viruses relative to virophages in one time step
x_virophage = 2*cell_radius
x_virus = sqrt(virophage_radius/virus_radius)*x_virophage
x_cell = sqrt(virophage_radius/cell_radius)*x_virophage
virus_burst_size = 100 #burst size for viruses
virophage_burst_size = 1000 #burst size for viruses

replication_cycle = 20 #number of time steps for a viral replication cycle
clearance_period = 100 #number of time steps for clearance of an exogenous virophage

programmed_cell_death = false #infected cells commit suicide imparing virus replication

# Initialise populations of agents
cell_number = 1 #number of cells in each organism {1,2,4,8,16}

## Initialisation

#function used in case of 16 cells
function find_beta(x,y,z,α,β,h)
    betas = Array{Float64}(undef,length(collect(0:0.0001:pi/2)),2)
    x1 = x+h*cos(α)*cos(β)
    y1 = y+h*sin(α)*cos(β)
    z1 = z+h*sin(β)
    x2 = x+h*cos(α)*cos(β+π/2)
    y2 = y+h*sin(α)*cos(β+π/2)
    z2 = z+h*sin(β+π/2)
    d1_2 = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
    n = 0
    for i in collect(0:0.0001:pi/2)
        n += 1
        new_β = β-i
        x11 = x+h*cos(α+π/2)*cos(new_β)
        y11 = y+h*sin(α+π/2)*cos(new_β)
        z11 = z+h*sin(new_β)
        d2_11 = sqrt((x2-x11)^2+(y2-y11)^2+(z2-z11)^2)
        betas[n,1] = new_β
        betas[n,2] = abs(d1_2-d2_11)
    end
    return betas[argmin(betas[:,2]),:][1]
end

#initialise cell population
cells = Array{cell,1}()
for i = 1:round(N0/cell_number)
    x = rand(Uniform(4*cell_radius,space-4*cell_radius))
    y = rand(Uniform(4*cell_radius,space-4*cell_radius))
    z= rand(Uniform(4*cell_radius,space-4*cell_radius))
    if cell_number == 1
        push!(cells,cell(1,i,[1,1],x,y,z))
    elseif cell_number == 2
        α = rand(Uniform(-2*pi,2*pi))
        β = rand(Uniform(-2*pi,2*pi))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α)*cos(β),y+cell_radius*sin(α)*cos(β),z+cell_radius*sin(β)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α)*cos(β+π),y+cell_radius*sin(α)*cos(β+π),z+cell_radius*sin(β+π)))
    elseif cell_number == 4
        α = rand(Uniform(-2*pi,2*pi))
        β = rand(Uniform(-2*pi,2*pi))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α)*cos(β),y+cell_radius*sin(α)*cos(β),z+cell_radius*sin(β)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α)*cos(β+π),y+cell_radius*sin(α)*cos(β+π),z+cell_radius*sin(β+π)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α)*cos(β+π/2),y+cell_radius*sin(α)*cos(β+π/2),z+cell_radius*sin(β+π/2)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α)*cos(β-π/2),y+cell_radius*sin(α)*cos(β-π/2),z+cell_radius*sin(β-π/2)))
    elseif cell_number == 8
        α = rand(Uniform(-2*pi,2*pi))
        β = rand(Uniform(-2*pi,2*pi))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α)*cos(β),y+cell_radius*sin(α)*cos(β),z+cell_radius*sin(β)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α)*cos(β-π/4),y+cell_radius*sin(α)*cos(β-π/4),z+cell_radius*sin(β-π/4)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α+π)*cos(-β),y+cell_radius*sin(α+π)*cos(-β),z+cell_radius*sin(-β)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α+π)*cos(-β+π/4),y+cell_radius*sin(α+π)*cos(-β+π/4),z+cell_radius*sin(-β+π/4)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α+π/2)*cos(β),y+cell_radius*sin(α+π/2)*cos(β),z+cell_radius*sin(β)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α+π/2)*cos(β-π/4),y+cell_radius*sin(α+π/2)*cos(β-π/4),z+cell_radius*sin(β-π/4)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α+3*π/2)*cos(-β),y+cell_radius*sin(α+3*π/2)*cos(-β),z+cell_radius*sin(-β)))
        push!(cells,cell(1,i,[1,1],x+cell_radius*cos(α+3*π/2)*cos(-β+π/4),y+cell_radius*sin(α+3*π/2)*cos(-β+π/4),z+cell_radius*sin(-β+π/4)))
    elseif cell_number == 16
        α = rand(Uniform(-2*pi,2*pi))
        β = rand(Uniform(-2*pi,2*pi))

        h = 2*cell_radius

        push!(cells,cell(1,i,[1,1],x+h*cos(α)*cos(β),y+h*sin(α)*cos(β),z+h*sin(β)))
        push!(cells,cell(1,i,[1,1],x+h*cos(α)*cos(β+π/2),y+h*sin(α)*cos(β+π/2),z+h*sin(β+π/2)))
        push!(cells,cell(1,i,[1,1],x+h*cos(α)*cos(β+π),y+h*sin(α)*cos(β+π),z+h*sin(β+π)))
        push!(cells,cell(1,i,[1,1],x+h*cos(α)*cos(β+3*π/2),y+h*sin(α)*cos(β+3*π/2),z+h*sin(β+3*π/2)))
        push!(cells,cell(1,i,[1,1],x+h/4*cos(α)*cos(β+π/2),y+h/4*sin(α)*cos(β+π/2),z+h/4*sin(β+π/2)))
        push!(cells,cell(1,i,[1,1],x+h/4*cos(α)*cos(β+3*π/2),y+h/4*sin(α)*cos(β+3*π/2),z+h/4*sin(β+3*π/2)))

        h2 = sqrt(h^2+(h*tan(pi/4))^2)

        push!(cells,cell(1,i,[1,1],x+h2*cos(α)*cos(β+π/4),y+h2*sin(α)*cos(β+π/4),z+h2*sin(β+π/4)))
        push!(cells,cell(1,i,[1,1],x+h2*cos(α)*cos(β+3*π/4),y+h2*sin(α)*cos(β+3*π/4),z+h2*sin(β+3*π/4)))
        push!(cells,cell(1,i,[1,1],x+h2*cos(α)*cos(β+5*π/4),y+h2*sin(α)*cos(β+5*π/4),z+h2*sin(β+5*π/4)))
        push!(cells,cell(1,i,[1,1],x+h2*cos(α)*cos(β+7*π/4),y+h2*sin(α)*cos(β+7*π/4),z+h2*sin(β+7*π/4)))

        new_β = find_beta(x,y,z,α,β,h)

        push!(cells,cell(1,i,[1,1],x+h*cos(α+π/2)*cos(new_β),y+h*sin(α+π/2)*cos(new_β),z+h*sin(new_β)))
        push!(cells,cell(1,i,[1,1],x+h*cos(α+π/2)*cos(new_β+π),y+h*sin(α+π/2)*cos(new_β+π),z+h*sin(new_β+π)))
        push!(cells,cell(1,i,[1,1],x+h2*cos(α+π/2)*cos(new_β+π/4),y+h2*sin(α+π/2)*cos(new_β+π/4),z+h2*sin(new_β+π/4)))
        push!(cells,cell(1,i,[1,1],x+h2*cos(α+π/2)*cos(new_β+3*π/4),y+h2*sin(α+π/2)*cos(new_β+3*π/4),z+h2*sin(new_β+3*π/4)))
        push!(cells,cell(1,i,[1,1],x+h2*cos(α+π/2)*cos(new_β+5*π/4),y+h2*sin(α+π/2)*cos(new_β+5*π/4),z+h2*sin(new_β+5*π/4)))
        push!(cells,cell(1,i,[1,1],x+h2*cos(α+π/2)*cos(new_β+7*π/4),y+h2*sin(α+π/2)*cos(new_β+7*π/4),z+h2*sin(new_β+7*π/4)))
    end
end

#initialise virus population
viruses = Array{virus,1}()
for i = 1:V0
    push!(viruses,virus(1,rand(Uniform(0,space)),rand(Uniform(0,space)),rand(Uniform(0,space))))
end

#initialise virophage population
virophages = Array{virophage,1}()
for i = 1:M0
    push!(virophages,virophage(1,rand(Uniform(0,space)),rand(Uniform(0,space)),rand(Uniform(0,space))))
end

## Functions

# Function definitions

function displacements(cell_array)
    #calculate number of unique organisms
    organisms = unique([cell_array[i].o for i in 1:length(cell_array)])
    #make dictionary with displacement vectors for each organism
    move = Dict{Int64,Tuple}()
    #get vectors displacement vectors for each organism
    for i in organisms
        #Gaussian random walk
        α = rand(Uniform(-2*pi,2*pi))
        β = rand(Uniform(-2*pi,2*pi))
        z = x_cell*sin(β)
        x = x_cell*cos(α)*cos(β)
        y = x_cell*sin(α)*cos(β)
        move[i] = (x,y,z)
    end
    return move
end

#check and ensure that displacements are bounded
function update_cell(cell,displacement_vectors)
    cell.x = cell.x + displacement_vectors[cell.o][1]
    cell.y = cell.y + displacement_vectors[cell.o][2]
    cell.z = cell.z + displacement_vectors[cell.o][3]
    #ensure that positions are bounded by space
    if cell.x < cell_radius
        cell.x = cell_radius
    elseif cell.x > space-cell_radius
        cell.x = space-cell_radius
    end
    if cell.y < cell_radius
        cell.y = cell_radius
    elseif cell.y > space-cell_radius
        cell.y = space-cell_radius
    end
    if cell.z < cell_radius
        cell.z = cell_radius
    elseif cell.z > space-cell_radius
        cell.z = space-cell_radius
    end
end

function update_virus(virus)
        α = rand(Uniform(-2*π,2*π))
        β = rand(Uniform(-2*π,2*π))
        virus.z = virus.z + x_virus*sin(β)
        virus.x = virus.x + x_virus*cos(α)*cos(β)
        virus.y = virus.y + x_virus*sin(α)*cos(β)
        #ensure that positions are bounded by space
        if virus.x < virus_radius
            virus.x = virus_radius
        elseif virus.x > space-virus_radius
            virus.x = space-virus_radius
        end
        if virus.y < virus_radius
            virus.y = virus_radius
        elseif virus.y > space-virus_radius
            virus.y = space-virus_radius
        end
        if virus.z < virus_radius
            virus.z = virus_radius
        elseif virus.z > space-virus_radius
            virus.z = space-virus_radius
    end
end

function update_virophage(virophage)
        α = rand(Uniform(-2*π,2*π))
        β = rand(Uniform(-2*π,2*π))
        virophage.z = virophage.z + x_virophage*sin(β)
        virophage.x = virophage.x + x_virophage*cos(α)*cos(β)
        virophage.y = virophage.y + x_virophage*sin(α)*cos(β)
        #ensure that positions are bounded by space
        if virophage.x < virophage_radius
            virophage.x = virophage_radius
        elseif virophage.x > space-virophage_radius
            virophage.x = space-virophage_radius
        end
        if virophage.y < virophage_radius
            virophage.y = virophage_radius
        elseif virophage.y > space-virophage_radius
            virophage.y = space-virophage_radius
        end
        if virophage.z < virophage_radius
            virophage.z = virophage_radius
        elseif virophage.z > space-virophage_radius
            virophage.z = space-virophage_radius
    end
end

function infect_cells_giant_virus(cell_array,virus_array,time)
    # get position vectors of cells
    cells_x = [cell_array[i].x for i = 1:length(cell_array)]
    cells_y = [cell_array[i].y for i = 1:length(cell_array)]
    cells_z = [cell_array[i].z for i = 1:length(cell_array)]
    #get position vectors of viruses
    viruses_x = [virus_array[i].x for i = 1:length(virus_array)]
    viruses_y = [virus_array[i].y for i = 1:length(virus_array)]
    viruses_z = [virus_array[i].z for i = 1:length(virus_array)]
    #postion arrays for NearestNeighbor Search
    cells_pos = vcat(cells_x',cells_y',cells_z')
    virus_pos = vcat(viruses_x',viruses_y',viruses_z')
    kdtree = KDTree(cells_pos)
    #find nearest cell neighbor to each virus, if distance < k, kill cell,
    #release virus
    ids, distances = knn(kdtree,virus_pos,1,true)
    kill_viruses = Array{Int64,1}()
    for i = 1:length(virus_array)
        if distances[i][1] < cell_radius+virus_radius
            push!(kill_viruses,i) #remove infecting virus
            if cell_array[ids[i][1]].s[1] == 1
                #cell infected by giant virus
                cell_array[ids[i][1]].s = [2,time]
            elseif cell_array[ids[i][1]].s[1] == 3
                #cell infected by giant virus and exogenous virophage
                cell_array[ids[i][1]].s = [5,time]
            elseif cell_array[ids[i][1]].s[1] == 4
                #cell infected by giant virus and endogenous virophage
                cell_array[ids[i][1]].s = [6,time]
            elseif cell_array[ids[i][1]].s[1] == 7
                #cell infected by endogenous, exogenous virophages
                #and giant virus
                cell_array[ids[i][1]].s = [8,time]
            end
        end
    end
    deleteat!(virus_array,unique(sort(kill_viruses)))
end

function infect_cells_virophages(cell_array,virophage_array,time)
    # get position vectors of cells
    cells_x = [cell_array[i].x for i = 1:length(cell_array)]
    cells_y = [cell_array[i].y for i = 1:length(cell_array)]
    cells_z = [cell_array[i].z for i = 1:length(cell_array)]
    #get position vectors of viruses
    virophages_x = [virophage_array[i].x for i = 1:length(virophage_array)]
    virophages_y = [virophage_array[i].y for i = 1:length(virophage_array)]
    virophages_z = [virophage_array[i].z for i = 1:length(virophage_array)]
    #postion arrays for NearestNeighbor Search
    cells_pos = vcat(cells_x',cells_y',cells_z')
    virophages_pos = vcat(virophages_x',virophages_y',virophages_z')
    kdtree = KDTree(cells_pos)
    #find nearest cell neighbor to each virus, if distance < k, kill cell,
    #release virus
    ids, distances = knn(kdtree,virophages_pos,1,true)
    kill_virophages = Array{Int64,1}()
    for i = 1:length(virophage_array)
        if distances[i][1] < cell_radius+virophage_radius
            push!(kill_virophages,i) #remove infecting virus
            if cell_array[ids[i][1]].s[1] == 1
                #cell infected by exogenous virophage
                cell_array[ids[i][1]].s = [3,time]
            elseif cell_array[ids[i][1]].s[1] == 2
                #cell infected by giant virus and exogenous virophage
                cell_array[ids[i][1]].s[1] = 5
            elseif cell_array[ids[i][1]].s[1] == 4
                #cell infected by endogenous, exogenous virophages
                cell_array[ids[i][1]].s[1] = 7
            elseif cell_array[ids[i][1]].s[1] == 6
                #cell infected by endogenous, exogenous virophage
                #and giant virus
                cell_array[ids[i][1]].s[1] = 8
            end
        end
    end
    deleteat!(virophage_array,unique(sort(kill_virophages)))
end

function integrate_virophages(cell_array,time)
    for my_cell in cell_array
        if my_cell.s[1] == 3
            yes_integrate = rand(Binomial(1,integration_p))
            if yes_integrate == 1
                my_cell.s = [4,time]
            end
        end
    end
end

function decay(array,time,k)
    remove = Array{Int64,1}()
    for i = 1:length(array)
        p_die = 1-exp(-k*(time-array[i].t))
        yes_die = rand(Binomial(1,p_die))
        if yes_die == 1
            push!(remove,i)
        end
    end
    deleteat!(array,unique(sort(remove)))
end

function cell_virus_virophage_dynamics(cell_array,virus_array,virophage_array,time)
    kill_cells = Array{Int64,1}()
    for i =1:length(cell_array)
        if cell_array[i].s[1] == 2 && time-cell_array[i].s[2] == replication_cycle
            #kill cell and release giant viruses
            x = cell_array[i].x
            y = cell_array[i].y
            z = cell_array[i].z
            push!(kill_cells,i)
            for i = 1:virus_burst_size
                push!(virus_array,virus(time,x,y,z))
            end
        elseif cell_array[i].s[1] == 3 && time-cell_array[i].s[2] == clearance_period
            #cell clears virophage
            cell_array[i].s = [1,time]
        elseif cell_array[i].s[1] == 5 && time-cell_array[i].s[2] == replication_cycle
            #produce virophages and giant viruses under inhibition
            x = cell_array[i].x
            y = cell_array[i].y
            z = cell_array[i].z
            push!(kill_cells,i)
            #print("inhibition!")
            for i = 1:round(virus_burst_size*f)
                push!(virus_array,virus(time,x,y,z))
            end
            for i = 1:virophage_burst_size
                push!(virophage_array,virophage(time,x,y,z))
            end
        elseif cell_array[i].s[1] == 6 && time-cell_array[i].s[2] == replication_cycle
            #produce virophages and giant viruses without inhibition
            #reactivation generates virophages but giant viruses are not inhibited
            x = cell_array[i].x
            y = cell_array[i].y
            z = cell_array[i].z
            push!(kill_cells,i)
            for i = 1:virus_burst_size
                push!(virus_array,virus(time,x,y,z))
            end
            for i = 1:virophage_burst_size
                push!(virophage_array,virophage(time,x,y,z))
            end
        elseif cell_array[i].s[1] == 7 && time-cell_array[i].s[2] == clearance_period
            #cell clears exogenous virophage
            cell_array[i].s = [1,time]
        elseif cell_array[i].s[1] == 8 && time-cell_array[i].s[2] == replication_cycle
            #produce virophages and giant viruses under inhibition
            x = cell_array[i].x
            y = cell_array[i].y
            z = cell_array[i].z
            push!(kill_cells,i)
            #print("inhibition!")
            for i = 1:round(virus_burst_size*f)
                push!(virus_array,virus(time,x,y,z))
            end
            for i = 1:virophage_burst_size
                push!(virophage_array,virophage(time,x,y,z))
            end
        end
    end
    deleteat!(cell_array,unique(sort(kill_cells)))
end

function pcd(cell_array)
    suicide = rand(Binomial(1,1/20)) #suicide probability in each time step 1/20
    #in a period of 20 time steps (virus replication cycle) probability is equal to 0.5
    kill_cells = Array{Int64,1}()
    for i =1:length(cell_array)
        if cell_array[i].s[1] == 2 && suicide == 1
            push!(kill_cells,i)
        elseif cell_array[i].s[1] == 5 && suicide == 1
            push!(kill_cells,i)
        elseif cell_array[i].s[1] == 6 && suicide == 1
            push!(kill_cells,i)
        end
    end
    deleteat!(cell_array,unique(sort(kill_cells)))
end

## Simulation

# Run simlutaion

t_max = 100000
df = Array{Int64}(undef,t_max,12)

#@gif 

for t = 1:t_max

    #println(t)

    if sum([1 for i = 1:length(cells) if cells[i].s[1] != 0]) > 1

    #update positions and number of cells
    displace = displacements(cells)
    for my_cell in cells
        update_cell(my_cell,displace)
    end
    #update position of viruses
    for my_virus in viruses
        update_virus(my_virus)
    end
    #update_position of virophages
    for my_virophage in virophages
        update_virophage(my_virophage)
    end
    #infection of cells by viruses
    if length(viruses) > 0
        infect_cells_giant_virus(cells,viruses,t)
    end
    #infection of cells by virophages
    if length(virophages) > 0
        infect_cells_virophages(cells,virophages,t)
    end
    #integration of virophage into the cell's genome
    integrate_virophages(cells,t)
    #virus decay
    decay(viruses,t,k1)
    #virophage decay
    decay(virophages,t,k2)
    #cell lysis, production of viruses and virophages
    cell_virus_virophage_dynamics(cells,viruses,virophages,t)
    #programmed cell-death
    if programmed_cell_death == true
        pcd(cells)
    end

    df[t,1] = t
    df[t,2] = sum([1 for i = 1:length(cells) if cells[i].s[1]==1])
    df[t,3] = sum([1 for i = 1:length(cells) if cells[i].s[1]==2])
    df[t,4] = sum([1 for i = 1:length(cells) if cells[i].s[1]==3])
    df[t,5] = sum([1 for i = 1:length(cells) if cells[i].s[1]==4])
    df[t,6] = sum([1 for i = 1:length(cells) if cells[i].s[1]==5])
    df[t,7] = sum([1 for i = 1:length(cells) if cells[i].s[1]==6])
    df[t,8] = sum([1 for i = 1:length(cells) if cells[i].s[1]==7])
    df[t,9] = sum([1 for i = 1:length(cells) if cells[i].s[1]==8])
    df[t,10] = sum([1 for i = 1:length(cells)]) #count cell population
    df[t,11] = length(viruses) #count virus population
    df[t,12] = length(virophages) #count virophage population

end

    #plot animation
    # get position vectors of cells
    cells_x = [cells[i].x for i = 1:length(cells)]
    cells_y = [cells[i].y for i = 1:length(cells)]
    cells_z = [cells[i].z for i = 1:length(cells)]
    #get position vectors of viruses
    viruses_x = [viruses[i].x for i = 1:length(viruses)]
    viruses_y = [viruses[i].y for i = 1:length(viruses)]
    viruses_z = [viruses[i].z for i = 1:length(viruses)]
    #get position vectors of virophages
    virophages_x = [virophages[i].x for i = 1:length(virophages)]
    virophages_y = [virophages[i].y for i = 1:length(virophages)]
    virophages_z = [virophages[i].z for i = 1:length(virophages)]
    #add frame to animated plot
    #=scatter(cells_x, cells_y, cells_z,
        xlim = (0, space),ylim = (0, space),zlim = (0, space),
        seriestype=:scatter, markersize = 3,leg=false, markercolor = :blue,
        title="time = "*string(t))
    scatter!(viruses_x, viruses_y, viruses_z,
        xlim = (0, space),ylim = (0, space),zlim = (0, space),
        seriestype=:scatter, markersize = 2,leg=false, markercolor = :red,
        title="time = "*string(t))
    scatter!(virophages_x, virophages_y, virophages_z,
        xlim = (0, space),ylim = (0, space),zlim = (0, space),
        seriestype=:scatter, markersize = 2,leg=false, markercolor = :green,
        title="time = "*string(t))=# 
end #every 10

results = vcat(["time" "uninfected" "infected_Gv" "infected_exo_Vp" "infected_endo_Vp" "infected_Gv_Vp" "infected_endo_Vp_Gv" "infected_endo_exo_Vp" "infected_Gv_endo_exo_Vp" "cells" "viruses" "virophages"],df)

writedlm("cvv_dynamics_cells1_neutral.csv",results,',')
