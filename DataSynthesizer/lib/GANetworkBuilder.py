from DataSynthesizer.lib.utils import mutual_information

import random
import copy

class GANetworkBuilder:

    def __init__(self, k, mutation_rate=None, source_genes=10, genepool_size=100, epochs=100, seed=0, sensi=None, target=None):
        self.k = k
        self.mutation_rate = mutation_rate
        self.source_genes = source_genes
        self.genepool_size = genepool_size
        self.epochs = epochs
        self.seed = seed
        self.sensi = sensi
        self.target = target

        self.fitness_map = {}

    def create_random_gene(self, df, n_cols):
        genome = []

        line0 = []
        for c in range(n_cols):
            line0.append(c)
        random.shuffle(line0)

        if self.sensi:
            temp_0=line0[0]
            target_ind = df.columns.tolist().index(self.target)
            line0[line0.index(target_ind)]=temp_0
            line0[0]=target_ind
            temp_1=line0[1]
            sensi_ind = df.columns.tolist().index(self.sensi)
            line0[line0.index(sensi_ind)]=temp_1
            line0[1]=sensi_ind

        genome.append(line0)

        for l in range(1,n_cols+1):
            clist=[x for x in line0 if line0.index(x)<line0.index(l-1)]
            if len(clist)>self.k:
                line = random.sample(clist,self.k)
            else:
                line = clist
                dlist=[x for x in line0 if x not in clist]
                line=line + random.sample(dlist,self.k-len(clist))
            genome.append(line)

        return genome


    def create_random_genepool(self, df, n):
        ret = []

        for i in range(self.genepool_size):
            ret.append(self.create_random_gene(df, len(df.columns)))

        return ret


    def get_mi(self, c, p, df):
        if self.fitness_map[c][p] == -1:
            self.fitness_map[c][p] = mutual_information(df[c], df[[p]])

        return self.fitness_map[c][p]


    def calc_fitness(self, genome, df):
        fitness = 0

        network = self.convert_to_network(genome, df)
        for (c, ps) in network:
            for p in ps:
                fitness += self.get_mi(c, p, df)
        return fitness


    def eval_genepool(self, genepool, df):
        evaluated_genepool = []

        for genome in genepool:
            fitness = self.calc_fitness(genome, df)
            evaluated_genepool.append((fitness, genome))

        evaluated_genepool.sort(reverse = True)
        return evaluated_genepool


    def crossover(self, genome_1, genome_2):
        n = random.randint(0, len(genome_1)-1)

        return genome_1[:n] + genome_2[n:]


    def flip(self, genome):

        if not self.mutation_rate: mr = 1/(len(genome)-1)
        else: mr = self.mutation_rate

        if not self.sensi:
                     for c1 in range(len(genome[0])):
                        if random.random() < mr:
                                c2=random.randint(0,len(genome[0])-1)
                                seq = genome[0][c1]
                                genome[0][c1] = genome[0][c2]
                                genome[0][c2] = seq
        else:
                     for c1 in range(2,len(genome[0])):
                        if random.random() < mr:
                                c2=random.randint(2,len(genome[0])-1)
                                seq = genome[0][c1]
                                genome[0][c1] = genome[0][c2]
                                genome[0][c2] = seq

        return genome


    def swap(self, genome):

        if not self.mutation_rate: mr = 1/(len(genome)-1)
        else: mr = self.mutation_rate


        for l in range(1, len(genome)):

                    glist=[x for x in genome[l] if genome[0].index(x)>genome[0].index(l-1)]
                    for ele in glist:
                        clist=[x for x in genome[0] if x not in genome[l] and genome[0].index(x)<genome[0].index(l-1)]
                        if len(clist)>0:
                            ind=genome[l].index(ele)
                            genome[l][ind]=random.choice(clist)

                    for c in range(0,len(genome[l])):
                        if random.random() < mr:
                            clist=[x for x in genome[0] if x not in genome[l] and genome[0].index(x)<genome[0].index(l-1)]
                            if len(clist)>0:
                                genome[l][c]=random.choice(clist)

        return genome


    def mutate(self, evaluated_genepool):

        evaluated_genepool = evaluated_genepool[:self.source_genes]
        new_genepool = []

        for eg in evaluated_genepool:
            new_genepool.append(eg[1])

        while len(new_genepool) < self.genepool_size:
            genome = copy.deepcopy(evaluated_genepool[random.randint(0, len(evaluated_genepool)-1)][1])

            if not self.mutation_rate:
                mr = 1/(len(genome)-1)
            else:
                mr = self.mutation_rate

            if random.random() < mr:
                genome_2 = copy.deepcopy(evaluated_genepool[random.randint(0, len(evaluated_genepool)-1)][1])
                genome = self.crossover(genome, genome_2)

            genome = self.flip(genome)
            genome = self.swap(genome)

            new_genepool.append(genome)

        return new_genepool


    def convert_to_network(self, genome, df):

        network = []

        if not self.sensi:
            for i, e in enumerate(genome[0]):
                if i == 0:
                    network.append((df.columns[e], []))
                elif i <= self.k:
                    network.append((df.columns[e], network[i-1][1]+[network[i-1][0]]))
                else:
                    prev_elements = genome[0][:i]
                    network.append((df.columns[e], []))
                    for e2 in genome[e+1]:
                        if e2 in prev_elements:
                            network[i][1].append(df.columns[e2])

            return network

        else:
            sensi_ind = df.columns.tolist().index(self.sensi)
            for i, e in enumerate(genome[0]):
                if i == 0:
                    network.append((df.columns[e], []))
                elif i <= self.k:
                    if i == 1 or i == 2:
                        network.append((df.columns[e],[self.target]))
                    else:
                        network.append((df.columns[e], network[i-1][1]+[network[i-1][0]]))
                else:
                    prev_elements = genome[0][:i]
                    network.append((df.columns[e], []))

                    for e2 in genome[e+1]:
                        if e2 in prev_elements and e2 != sensi_ind:
                            network[i][1].append(df.columns[e2])

            return network



    def ga_network(self, df):
        str_df = df.astype(str, copy=False)
        genepool = self.create_random_genepool(str_df, self.genepool_size)

        random.seed(self.seed)

        self.fitness_map = {}
        for c in str_df.columns:
            self.fitness_map[c]={}
            for p in str_df.columns:
                self.fitness_map[c][p] = -1

        for epoch in range(0,self.epochs):
            evaluated_genepool = self.eval_genepool(genepool, str_df)
            print("Epoch {}/{}".format(epoch+1, self.epochs))
            #print(evaluated_genepool[0][0])
            print(evaluated_genepool[0])
            print("-----\n")
            genepool = self.mutate(evaluated_genepool)

        return self.convert_to_network(evaluated_genepool[0][1], str_df)
