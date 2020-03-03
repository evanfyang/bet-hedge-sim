# Stochastic switching population simulation
# Plot the mean value of the switching rate

# Selection is not working, find out what is going on
# -> Investigate fitness values, figure out why it is not changing 
# When environment is switching, members with phenotype matching the environment should be more common
# -> Investigate how fitness values are being calculated 

using Plots
using Distributions
import Base.copy, Base.copy!

#pyplot()

#using Pkg
#Pkg.add("PyPlot")


# Individuals of a Population.
mutable struct Individual
    phenotype::Int64            # 0 or 1
    genotype::Float64           # switching rate
    fitness::Float64            # survival

    # Constructor
    # Individual() = new(Int64, Float64, Float64)
    # Individual(phenotype) = new(phenotype, Float64, Float64)
    # Individual(phenotype, genotype) = new(phenotype, genotype, Float64)
    Individual(phenotype, genotype, fitness) = new(phenotype, genotype, fitness)
end

# Population - collection of Individuals with population-level properties.
mutable struct Population
    size::Int64
    env_state::Int64
    env_switch_rate::Int64
    init_switch_rate::Float64
    selec_coeff_0::Float64
    selec_coeff_1::Float64
    members::Array{Individual,1}
    members_prev::Array{Individual,1}
    fitness::Array{Float64,1}

    # Constructor
    function Population(size::Int64, env_state::Int64, env_switch_rate::Int64,
                        init_switch_rate::Float64, selec_coeff_0::Float64,
                        selec_coeff_1::Float64)
        members = Array{Individual,1}(undef, size)
        # Randomly set phenotype to 0 or 1 for each
        # individual in initial population.
        for i=1:size
            p = rand(Uniform(0, 1))
            if p <= 0.5
                p = 0
            else
                p = 1
            end
            ind = Individual(p, init_switch_rate, 1.0)
            members[i] = ind
        end
        members_prev = copy(members)
        new(size, env_state, env_switch_rate, init_switch_rate, selec_coeff_0,
        selec_coeff_1, members, members_prev, ones(size))
    end
end

# copy inplace method for Individual
function copy!(i::Individual, j::Individual)
    i.genotype = j.genotype
    i.phenotype = j.phenotype
    i.fitness = j.fitness
end

# copy inplace method for Array of Individuals
function copy!(x::Array{Individual,1}, y::Array{Individual,1})
    if (length(x) != length(y))
        error("Arrays must have equal length")
    end

    for i=1:length(x)
        copy!(x[i], y[i])
    end
end

# function BD(μ, σ)
#     a = μ * (μ - μ^2 - σ^2) / σ^2
#     b = (1 - μ) * (μ - μ^2 - σ^2) / σ^2
#     @assert a > zero(a) && b > zero(b) "BD: $μ, $σ, $a, $b"
#     return Beta(a,b)
# end

# function mutate(parent_rate, mutation_mult)
#     parent_rate >= 1.0 - eps(Float64) ? parent_rate = 1.0 - 2*eps(Float64) : nothing
#     parent_rate <= eps(Float64) ? parent_rate = 2*eps(Float64) : nothing
#     @assert !isnan(parent_rate) && !isnan(mutation_mult) "mutate: $parent_rate, $mutation_mult"
#     beta_std = mutation_mult * sqrt(parent_rate * (1-parent_rate))
#     return rand(BD(parent_rate, beta_std))
# end

function mutate(parent, sigma)
    @assert parent > 0 && parent < 1 "parent value $parent out of bounds"
    if rand() < 0.5
        return parent - rand(Truncated(Exponential(sigma), 0, parent))
    else
        return parent + rand(Truncated(Exponential(sigma), 0, 1-parent))
    end
end

# Main life cycle function
function next_gen(pop::Population, mutation_rate::Float64, switch_env::Bool)
    # Save previous population
    copy!(pop.members_prev, pop.members)

    # Update fitness of each individual
    sum_fit = 0;
    for i=1:pop.size
        if pop.members[i].phenotype == pop.env_state
            pop.members[i].fitness = 1
        else
            if pop.members[i].phenotype == 0
                pop.members[i].fitness = 1 - pop.selec_coeff_1
            else
                pop.members[i].fitness = 1 - pop.selec_coeff_0
            end
        end
        pop.fitness[i] = pop.members[i].fitness
        sum_fit += pop.members[i].fitness
    end
    mean_fit = sum_fit / pop.size
    # Normalize fitness for each individual
    for i=1:pop.size
        pop.members[i].fitness = pop.members[i].fitness / mean_fit
        pop.fitness[i] = pop.fitness[i] / sum_fit
    end

    # Survival and Reproduction
    fit_dist = Categorical(pop.fitness)
    for i=1:pop.size
        parent = rand(fit_dist) # should not be uniform, should sample based on fitness
        copy!(pop.members[i], pop.members_prev[parent])
        # Do phenotype switching based on parent's switching
        # do not sample from mean fitness, sample from categorical distribution
        offspring_switch_phenotype = rand(Uniform(0, 1))
        if offspring_switch_phenotype <= pop.members_prev[parent].genotype
            if pop.members[i].phenotype == 0
                pop.members[i].phenotype = 1
            else
                pop.members[i].phenotype = 0
            end
        end
        # Mutate switching rate using mutation probability
        # rand(Truncated(Normal(mu, sigma), 0, 1))
        # sigma should be chosen to be relatively small, say 0.01 or 0.05
        # mu is the parental switch rate.
        offspring_mutate_switch_rate = rand(Uniform(0, 1))
        if offspring_mutate_switch_rate <= mutation_rate
            ## start with 0.05, run 100 times, see last value, average it.
            pop.members[i].genotype = mutate(pop.members_prev[parent].genotype, 0.01)
        end
    end

    # Update environmental state
    if switch_env
        if pop.env_state == 0
            pop.env_state = 1
        else
            pop.env_state = 0
        end
    end
end

# Main function
function run_sim(num_generations::Int64, popSize::Int64, env_switch_rate::Int64,
                 init_switch_rate::Float64, mutation_rate::Float64,
                 selec_coeff_1::Float64, selec_coeff_2::Float64)
    pop = Population(popSize, 0, env_switch_rate, init_switch_rate, selec_coeff_1, selec_coeff_2)

    gen_mean_genotype = Array{Float64,1}(undef, num_generations)
    env_states = zeros(num_generations)
    fitness = zeros(num_generations, popSize)
    genotypes = zeros(num_generations, popSize)
    phenotypes = zeros(num_generations, popSize)

    # Switch the environmental state every num_generations
    for i=1:num_generations
        if i % pop.env_switch_rate == 0
            switch_env = true
        else
            switch_env = false
        end
        env_states[i] = pop.env_state
        next_gen(pop, mutation_rate, switch_env)
        # Calculate mean switching rate for each generation
        sum_switch = 0;
        for j=1:pop.size
            sum_switch += pop.members[j].genotype
        end
        mean_switch = sum_switch / pop.size
        gen_mean_genotype[i] = mean_switch
        fitness[i,:] = pop.fitness
        genotypes[i,:] = [pop.members[i].genotype for i in 1:pop.size]
        phenotypes[i,:] = [pop.members[i].phenotype for i in 1:pop.size]
    end

    heatmap(pop.fitness, pop.fitness)

    return gen_mean_genotype, fitness, env_states, genotypes, phenotypes
end

function main(args)
    #run_sim(100, 1000, 5, 0.5, 0.7, 0.4, 0.4)
    #TODO: find a way to pass arguments to run_sim from command line.
end

## switch every 20 generations, should get 0.5
## To run, cmd a and then cmd enter or ctl enter
main()
