# Stochastic switching population simulation
# Plot the mean value of the switching rate
# 


using Distributions

# Individuals of a Population.
struct Individual
    age::Int64
    fitness::Float64            # survival
    phenotype::Int64            # 0 or 1
    genotype::Float64           # switching rate 

    # Constructor
    Individual() = new(0.0, Float64, Int64, Float64)
    Individual(age) = new(age, Float64, Int64, Float64)
    Individual(age, phenotype) = new(age, phenotype, Int64, Float64)
    Individual(age, phenotype, genotype) = new(age, phenotype, genotype, Float64)
    Individual(age, phenotype, genotype, fitness) = new(age, phenotype, genotype, fitness)

end

# Population - collection of Individuals with population-level properties.
struct Population
    size::Int64
    env_state::Int64
    env_switch_rate::Int64
    init_switch_rate::Float64
    members::Array{Individual,1}
    members_prev::Array{Individual,1}

    # Constructor
    function Population(size::Int64, env_state::Int64, env_switch_rate::Int64, 
                        init_switch_rate::Float64)
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
            ind = Individual(0, p, init_switch_rate)
            members[i] = ind
        end
        members_prev = copy(members)
    end
end

function mutate_switch_rate(offspring::Individual, parent::Individual)
    # Unsure of how to mutate switching rate...
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

            pop.members[i].fitness = 1 - # should be selection coefficient, not pop.members[i].genotype
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
            if offspring_switch_phenotype <= pop.members_prev[parent]).genotype
                if pop.members[i].phenotype == 0
                    pop.members[i].phenotype = 1
                else
                    pop.members[i].phenotype = 0
                end
            end
            # Mutate switching rate using mutation probability
            offspring_mutate_switch_rate = rand(Uniform(0, 1))
            if offspring_mutate_switch_rate <= mutation_rate
                mutate_switch_rate(pop.members[i], pop.members_prev[parent]) 
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
# need to pass in two selection coefficient
function run_sim(num_generations::Int64, popSize::Int64, env_switch_rate::Float64, 
                 init_switch_rate::Float64, mutation_rate::Float64)
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