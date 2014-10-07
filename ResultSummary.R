#Combine Result Comparision######


Megan_ni_4m11_500mnds_1nrep_100_out <- read.csv("Megan_ni_4m11_500mnds_1nrep_100_out.csv")
Megan_ni_4m11_500mnds_2nrep_100_out <- read.csv("Megan_ni_4m11_500mnds_2nrep_100_out.csv")
Megan_ni_4m11_1000mnds_1nrep_100_out <- read.csv("Megan_ni_4m11_1000mnds_1nrep_100_out.csv")
Megan_ni_4m11_1000mnds_2nrep_100_out <- read.csv("Megan_ni_4m11_1000mnds_2nrep_100_out.csv")
Megan_ni_4m11_2000mnds_1nrep_100_out <- read.csv("Megan_ni_4m11_2000mnds_1nrep_100_out.csv")
Megan_ni_4m11_2000mnds_2nrep_100_out <- read.csv("Megan_ni_4m11_2000mnds_2nrep_100_out.csv")
Megan_ni_4m11_3000mnds_1nrep_100_out <- read.csv("Megan_ni_4m11_3000mnds_1nrep_100_out.csv")
Megan_ni_4m11_3000mnds_2nrep_100_out <- read.csv("Megan_ni_4m11_3000mnds_2nrep_100_out.csv")


Megan_ni_10m11_500mnds_1nrep_100_out <- read.csv("Megan_ni_10m11_500mnds_1nrep_100_out.csv")
Megan_ni_10m11_500mnds_2nrep_100_out <- read.csv("Megan_ni_10m11_500mnds_2nrep_100_out.csv")
Megan_ni_10m11_1000mnds_1nrep_100_out <- read.csv("Megan_ni_10m11_1000mnds_1nrep_100_out.csv")

Megan_out <- rbind(Megan_ni_4m11_500mnds_1nrep_100_out, 
             Megan_ni_4m11_500mnds_2nrep_100_out,
             Megan_ni_4m11_1000mnds_1nrep_100_out,
             Megan_ni_4m11_1000mnds_2nrep_100_out,
             Megan_ni_4m11_2000mnds_1nrep_100_out,
             Megan_ni_4m11_2000mnds_2nrep_100_out,
             Megan_ni_4m11_3000mnds_1nrep_100_out,
             Megan_ni_4m11_3000mnds_2nrep_100_out,
             Megan_ni_10m11_500mnds_1nrep_100_out,
             Megan_ni_10m11_500mnds_2nrep_100_out,
             Megan_ni_10m11_1000mnds_1nrep_100_out)

row.names(Megan_out) <- c("ni=4, m11= 500, mnds = 1", 
                          "1",
                         "ni=4, m11= 500, mnds = 2", 
                         "2",
                         "ni=4, m11= 1000, mnds = 1", 
                         "3",
                         "ni=4, m11= 1000, mnds = 2",
                         "4",
                         "ni=4, m11= 2000, mnds = 1", 
                         "5",
                         "ni=4, m11= 2000, mnds = 2", 
                         "6",
                         "ni=4, m11= 3000, mnds = 1", 
                         "7",
                         "ni=4, m11= 3000, mnds = 2", 
                         "8",
                         "ni=10, m11= 500, mnds = 1", 
                         "9",
                         "ni=10, m11= 500, mnds = 2", 
                         "10",
                         "ni=10, m11= 1000, mnds = 1",
                         "11")
write.csv(Megan_out, "Megan_out.csv", row.names = T)


phillips_n_2000pi0_0.95mua_1nrep_100_out <- read.csv("phillips_n_2000pi0_0.95mua_1nrep_100_out.csv")
phillips_n_2000pi0_0.95mua_2nrep_100_out <- read.csv("phillips_n_2000pi0_0.95mua_2nrep_100_out.csv")
phillips_n_2000pi0_0.95mua_3nrep_100_out <- read.csv("phillips_n_2000pi0_0.95mua_3nrep_100_out.csv")
phillips_n_2000pi0_0.95mua_4nrep_100_out <- read.csv("phillips_n_2000pi0_0.95mua_4nrep_100_out.csv")

phillips_n_2000pi0_0.9mua_1nrep_100_out <- read.csv("phillips_n_2000pi0_0.9mua_1nrep_100_out.csv")
phillips_n_2000pi0_0.9mua_2nrep_100_out <- read.csv("phillips_n_2000pi0_0.9mua_2nrep_100_out.csv")
phillips_n_2000pi0_0.9mua_3nrep_100_out <- read.csv("phillips_n_2000pi0_0.9mua_3nrep_100_out.csv")
phillips_n_2000pi0_0.9mua_4nrep_100_out <- read.csv("phillips_n_2000pi0_0.9mua_4nrep_100_out.csv")

phillips_n_2000pi0_0.8mua_1nrep_100_out <- read.csv("phillips_n_2000pi0_0.8mua_1nrep_100_out.csv")
phillips_n_2000pi0_0.8mua_2nrep_100_out <- read.csv("phillips_n_2000pi0_0.8mua_2nrep_100_out.csv")
phillips_n_2000pi0_0.8mua_3nrep_100_out <- read.csv("phillips_n_2000pi0_0.8mua_3nrep_100_out.csv")
phillips_n_2000pi0_0.8mua_4nrep_100_out <- read.csv("phillips_n_2000pi0_0.8mua_4nrep_100_out.csv")


phillips_n_2000pi0_0.7mua_1nrep_100_out <- read.csv("phillips_n_2000pi0_0.7mua_1nrep_100_out.csv")
phillips_n_2000pi0_0.7mua_2nrep_100_out <- read.csv("phillips_n_2000pi0_0.7mua_2nrep_100_out.csv")
phillips_n_2000pi0_0.7mua_3nrep_100_out <- read.csv("phillips_n_2000pi0_0.7mua_3nrep_100_out.csv")
phillips_n_2000pi0_0.7mua_4nrep_100_out <- read.csv("phillips_n_2000pi0_0.7mua_4nrep_100_out.csv")


phillips_n_5000pi0_0.95mua_1nrep_100_out <- read.csv("phillips_n_5000pi0_0.95mua_1nrep_100_out.csv")
phillips_n_5000pi0_0.95mua_2nrep_100_out <- read.csv("phillips_n_5000pi0_0.95mua_2nrep_100_out.csv")
phillips_n_5000pi0_0.95mua_3nrep_100_out <- read.csv("phillips_n_5000pi0_0.95mua_3nrep_100_out.csv")
phillips_n_5000pi0_0.95mua_4nrep_100_out <- read.csv("phillips_n_5000pi0_0.95mua_4nrep_100_out.csv")

phillips_n_5000pi0_0.9mua_1nrep_100_out <- read.csv("phillips_n_5000pi0_0.9mua_1nrep_100_out.csv")
phillips_n_5000pi0_0.9mua_2nrep_100_out <- read.csv("phillips_n_5000pi0_0.9mua_2nrep_100_out.csv")
phillips_n_5000pi0_0.9mua_3nrep_100_out <- read.csv("phillips_n_5000pi0_0.9mua_3nrep_100_out.csv")
phillips_n_5000pi0_0.9mua_4nrep_100_out <- read.csv("phillips_n_5000pi0_0.9mua_4nrep_100_out.csv")

phillips_n_5000pi0_0.8mua_1nrep_100_out <- read.csv("phillips_n_5000pi0_0.8mua_1nrep_100_out.csv")
phillips_n_5000pi0_0.8mua_2nrep_100_out <- read.csv("phillips_n_5000pi0_0.8mua_2nrep_100_out.csv")
phillips_n_5000pi0_0.8mua_3nrep_100_out <- read.csv("phillips_n_5000pi0_0.8mua_3nrep_100_out.csv")
phillips_n_5000pi0_0.8mua_4nrep_100_out <- read.csv("phillips_n_5000pi0_0.8mua_4nrep_100_out.csv")

phillips_n_5000pi0_0.7mua_4nrep_100_out
phillips_n_5000pi0_0.7mua_1nrep_100_out <- read.csv("phillips_n_5000pi0_0.7mua_1nrep_100_out.csv")
phillips_n_5000pi0_0.7mua_2nrep_100_out <- read.csv("phillips_n_5000pi0_0.7mua_2nrep_100_out.csv")
phillips_n_5000pi0_0.7mua_3nrep_100_out <- read.csv("phillips_n_5000pi0_0.7mua_3nrep_100_out.csv")
#phillips_n_5000pi0_0.7mua_4nrep_100_out <- read.csv("phillips_n_5000pi0_0.7mua_4nrep_100_out.csv")

####
phillips_n_10000pi0_0.95mua_1nrep_100_out <- read.csv("phillips_n_10000pi0_0.95mua_1nrep_100_out.csv")
phillips_n_10000pi0_0.95mua_2nrep_100_out <- read.csv("phillips_n_10000pi0_0.95mua_2nrep_100_out.csv")
phillips_n_10000pi0_0.95mua_3nrep_100_out <- read.csv("phillips_n_10000pi0_0.95mua_3nrep_100_out.csv")
phillips_n_10000pi0_0.95mua_4nrep_100_out <- read.csv("phillips_n_10000pi0_0.95mua_4nrep_100_out.csv")

phillips_n_10000pi0_0.9mua_1nrep_100_out <- read.csv("phillips_n_10000pi0_0.9mua_1nrep_100_out.csv")
phillips_n_10000pi0_0.9mua_2nrep_100_out <- read.csv("phillips_n_10000pi0_0.9mua_2nrep_100_out.csv")
phillips_n_10000pi0_0.9mua_3nrep_100_out <- read.csv("phillips_n_10000pi0_0.9mua_3nrep_100_out.csv")
phillips_n_10000pi0_0.9mua_4nrep_100_out <- read.csv("phillips_n_10000pi0_0.9mua_4nrep_100_out.csv")

phillips_n_10000pi0_0.8mua_1nrep_100_out <- read.csv("phillips_n_10000pi0_0.8mua_1nrep_100_out.csv")
phillips_n_10000pi0_0.8mua_2nrep_100_out <- read.csv("phillips_n_10000pi0_0.8mua_2nrep_100_out.csv")
phillips_n_10000pi0_0.8mua_3nrep_100_out <- read.csv("phillips_n_10000pi0_0.8mua_3nrep_100_out.csv")
phillips_n_10000pi0_0.8mua_4nrep_100_out <- read.csv("phillips_n_10000pi0_0.8mua_4nrep_100_out.csv")


phillips_n_10000pi0_0.7mua_1nrep_100_out <- read.csv("phillips_n_10000pi0_0.7mua_1nrep_100_out.csv")
phillips_n_10000pi0_0.7mua_2nrep_100_out <- read.csv("phillips_n_10000pi0_0.7mua_2nrep_100_out.csv")
phillips_n_10000pi0_0.7mua_3nrep_100_out <- read.csv("phillips_n_10000pi0_0.7mua_3nrep_100_out.csv")



Daisy_out <- rbind(phillips_n_2000pi0_0.95mua_1nrep_100_out,
                   phillips_n_2000pi0_0.95mua_2nrep_100_out,
                   phillips_n_2000pi0_0.95mua_3nrep_100_out,
                   phillips_n_2000pi0_0.95mua_4nrep_100_out,
                   
                   phillips_n_2000pi0_0.9mua_1nrep_100_out,
                   phillips_n_2000pi0_0.9mua_2nrep_100_out,
                   phillips_n_2000pi0_0.9mua_3nrep_100_out,
                   phillips_n_2000pi0_0.9mua_4nrep_100_out,
                   
                   phillips_n_2000pi0_0.8mua_1nrep_100_out,
                   phillips_n_2000pi0_0.8mua_2nrep_100_out,
                   phillips_n_2000pi0_0.8mua_3nrep_100_out,
                   phillips_n_2000pi0_0.8mua_4nrep_100_out,
                   
                   phillips_n_2000pi0_0.7mua_1nrep_100_out,
                   phillips_n_2000pi0_0.7mua_2nrep_100_out,
                   phillips_n_2000pi0_0.7mua_3nrep_100_out,
                   phillips_n_2000pi0_0.7mua_4nrep_100_out,
                   
                   phillips_n_5000pi0_0.95mua_1nrep_100_out,
                   phillips_n_5000pi0_0.95mua_2nrep_100_out,
                   phillips_n_5000pi0_0.95mua_3nrep_100_out,
                   phillips_n_5000pi0_0.95mua_4nrep_100_out,
                   
                   phillips_n_5000pi0_0.9mua_1nrep_100_out,
                   phillips_n_5000pi0_0.9mua_2nrep_100_out,
                   phillips_n_5000pi0_0.9mua_3nrep_100_out,
                   phillips_n_5000pi0_0.9mua_4nrep_100_out,
                   
                   phillips_n_5000pi0_0.8mua_1nrep_100_out,
                   phillips_n_5000pi0_0.8mua_2nrep_100_out,
                   phillips_n_5000pi0_0.8mua_3nrep_100_out,
                   phillips_n_5000pi0_0.8mua_4nrep_100_out,
                   
                   phillips_n_5000pi0_0.7mua_1nrep_100_out,
                   phillips_n_5000pi0_0.7mua_2nrep_100_out,
                   phillips_n_5000pi0_0.7mua_3nrep_100_out,
                   #phillips_n_5000pi0_0.7mua_4nrep_100_out,
                   
                   phillips_n_10000pi0_0.95mua_1nrep_100_out,
                   phillips_n_10000pi0_0.95mua_2nrep_100_out,
                   phillips_n_10000pi0_0.95mua_3nrep_100_out,
                   phillips_n_10000pi0_0.95mua_4nrep_100_out,
                   
                   phillips_n_10000pi0_0.9mua_1nrep_100_out,
                   phillips_n_10000pi0_0.9mua_2nrep_100_out,
                   phillips_n_10000pi0_0.9mua_3nrep_100_out,
                   phillips_n_10000pi0_0.9mua_4nrep_100_out,
                   
                   phillips_n_10000pi0_0.8mua_1nrep_100_out,
                   phillips_n_10000pi0_0.8mua_2nrep_100_out,
                   phillips_n_10000pi0_0.8mua_3nrep_100_out,
                   phillips_n_10000pi0_0.8mua_4nrep_100_out,
                   
                   phillips_n_10000pi0_0.7mua_1nrep_100_out,
                   phillips_n_10000pi0_0.7mua_2nrep_100_out,
                   phillips_n_10000pi0_0.7mua_3nrep_100_out
                   )

write.csv(Daisy_out, file = "Daisy_out.csv")

# 2000, p0 4 mua 4
# 5000
# 10000
