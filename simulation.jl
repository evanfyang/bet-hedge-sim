# Stochastic switching population simulation
# Plot the mean value of the switching rate

using Distributions
import Base.copy, Base.copy!

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
        selec_coeff_1, members, members_prev)
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

function BD(μ, σ)
    a = μ * (μ - μ^2 - σ^2) / σ^2
    b = (1 - μ) * (μ - μ^2 - σ^2) / σ^2
    return Beta(a,b)
end

function mutate(parent_rate, mutation_mult)
    beta_std = mutation_mult * sqrt(parent_rate * (1-parent_rate))
    return rand(BD(parent_rate, beta_std))
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
        sum_fit += pop.members[i].fitness
    end
    mean_fit = sum_fit / pop.size
    # Normalize fitness for each individual
    for i=1:pop.size
        pop.members[i].fitness = pop.members[i].fitness / mean_fit
    end
    
    # Survival and Reproduction
    for i=1:pop.size
        parent = rand(1:pop.size) # should not be uniform, should sample based on fitness 
        # Do phenotype switching based on parent's switching
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
    # Switch the environmental state every num_generations
    for i=1:num_generations
        if pop.env_switch_rate % i == 0
            switch_env = true
        else
            switch_env = false
        end
        next_gen(pop, mutation_rate, switch_env)
        print(pop.members)
    end
end

function main()
    run_sim(1000, 500, 10, 0.5, 0.7, 0.4, 0.3)
end

main()