include("anchor_topic_profile.jl")

# Load the data
(words,Q) = load_uci("enron"; min_tf=0, min_tfidf=0)

Profile.init(10^6, 0.1)

# Run the topic mining
mine_topics(Q,100)

@profile begin
mine_topics(Q,100)
end;

Profile.print(format=:flat)
