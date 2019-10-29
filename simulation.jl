# Stochastic switching population simulation
# Plot the mean value of the switching rate

using Distributions

# Individuals of a Population.
struct Individual
    fitness::Float64            # survival
    phenotype::Int64            # 0 or 1
    genotype::Float64           # switching rate 

    # Constructor
    Individual() = new(Float64, Int64, Float64)
    Individual(phenotype) = new(phenotype, Int64, Float64)
    Individual(phenotype, genotype) = new(phenotype, genotype, Float64)
    Individual(phenotype, genotype, fitness) = new(phenotype, genotype, fitness)
end

# Population - collection of Individuals with population-level properties.
struct Population
    env_state::Int64
    env_switch_rate::Int64
    init_switch_rate::Float64
    selec_coeff_1::Float64
    selec_coeff_2::Float64
    members::Array{Individual,1}
    members_prev::Array{Individual,1}

    # Constructor
    function Population(size::Int64, env_state::Int64, env_switch_rate::Int64, 
                        init_switch_rate::Float64, selec_coeff_1::Float64,
                        selec_coeff_2::Float64)
        members = Array{Individual}(size)
        # Randomly set phenotype to 0 or 1 for each 
        # individual in initial population.
        for i=1:size
            p = rand(Uniform(0, 1))
            if p <= 0.5
                p = 0
            else
                p = 1
            end
            ind = Individual(p, init_switch_rate, 1)
            members[i] = ind
        end
        members_prev = copy(members)
    end
end

# Main life cycle function
function next_gen(pop::Population, mutation_rate::Float64, switch_env::Bool)
    # Save previous population 
    pop.members_prev = copy(pop.members)
    
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
        sum += pop.members[i].fitness
    end
    mean_fit = sum_fit / pop.size
    # Normalize fitness for each individual
    for i=1:pop.size
        pop.members[i].fitness = pop.members[i].fitness / mean_fit
    
    # Survival and Reproduction
    for i=1:pop.size
        if pop.members[i].fitness > rand()
            # Individual survives and ages
            pop.members[i].age += 1 # get rid of age.
        else
            # Individual dies and is replaced by random new born 
            pop.members[i].age = 0
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
                pop.members[i].genotype = rand(Truncated(Normal(mu, 0.01), 0, 1))
            end
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
function run_sim(num_generations::Int64, popSize::Int64, env_switch_rate::Float64, 
                 init_switch_rate::Float64, mutation_rate::Float64, 
                 selec_coeff_1::Float64, selec_coeff_2::Float64)
    pop = Population(popSize, 0, env_switch_rate, init_switch_rate)
    # Switch the environmental state every num_generations
    for i=1:num_generations
        if pop.env_switch_rate % i == 0
            switch_env = true
        else
            switch_env = false
        end
        next_gen(pop, mutation_rate, switch_env)
    end
end
