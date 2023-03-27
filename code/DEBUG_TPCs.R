
data_df <- filter(data.in.TPC, mosquito_species == mosquito_in, trait.name == trait_in)
plot(data_df$T, data_df$trait)

# sample_index = sample(300,1)
for (sample_index in 1:dim(temp_sample)[1]) {
(c = filter(temp_sample, sample_num == sample_index)$c)
(T0 = filter(temp_sample, sample_num == sample_index)$T0)
(Tm = filter(temp_sample, sample_num == sample_index)$Tm)

test_func <- Quadratic(c,T0,Tm)



lines(seq(0,50, by = 0.1), test_func(seq(0,50, by = 0.1)))
}