# -- Driver for the anchor words
#using IProfile
#using Profile

include("anchor_topic_profile.jl")

# Load the data
words = readcsv("yelp/words.csv")
Q = readcsv("yelp/Q.csv")

# Run the topic mining
mine_topics(Q)

#Q = readcsv("yelp/Q.csv")

@profile begin
mine_topics(Q)
mine_topics(Q)
mine_topics(Q)
mine_topics(Q)
end;

Profile.print(format=:flat)

# Write the results files
#write_topics("topics.txt", words, p, r, A, TW)
