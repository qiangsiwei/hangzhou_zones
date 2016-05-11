import re
import json
import math
import random
import fileinput
import collections

class LDASampler(object):
    def __init__(self, docs=None, num_topics=None, alpha=0.1, beta=0.1, state=None):
        if state:
            self.__dict__ = json.loads(state)
        else:
            self.a, self.b, self.T = float(alpha), float(beta), int(num_topics)
            self.docs = docs
            # self.vocab = list(set(word for doc in self.docs for word in doc))
            self.vocab = [str(i)+ct for ct in [",1,E",",1,D",",0,E",",0,D"] for i in xrange(24)]
            self.D, self.W = len(self.docs), len(self.vocab)  
            self.to_int = {word: w for (w, word) in enumerate(self.vocab)}
            self.nt = [0] * self.T
            self.nd = [len(doc) for doc in self.docs]
            self.nwt = [[0] * self.T for _ in self.vocab]
            self.ndt = [[0] * self.T for _ in self.docs]
            self.assignments = []
            for d, doc in enumerate(docs):
                for i, word in enumerate(doc):
                    w = self.to_int[word]
                    t = random.randint(0, self.T - 1)
                    z = [d, i, w, t]
                    self.assignments.append(z)
            for z in self.assignments:
                d, _, w, t = z
                self.nt[t] += 1
                self.nwt[w][t] += 1
                self.ndt[d][t] += 1
    def wordmap(self):
        return self.to_int
    def next(self):
        for z in self.assignments:
            self.sample(z)
    def sample(self, z):
        d, _, w, old_t = z
        self.nt[old_t] -= 1
        self.ndt[d][old_t] -= 1
        self.nwt[w][old_t] -= 1
        unnorm_ps = [self.f(d, w, t) for t in range(self.T)]
        r = random.random() * sum(unnorm_ps)
        new_t = self.T - 1
        for i in range(self.T):
            r = r - unnorm_ps[i]
            if r < 0:
                new_t = i
                break
        z[3] = new_t
        self.nt[new_t] += 1
        self.ndt[d][new_t] += 1
        self.nwt[w][new_t] += 1
    def f(self, d, w, t):
        return ((self.ndt[d][t] + self.a) * (self.nwt[w][t] + self.b) / (self.nt[t] + self.W * self.b))
    def pw_z(self, w, t):
        return (self.nwt[w][t] + self.b) / (self.nt[t] + self.W * self.b)
    def pz_d(self, d, t):
        return (self.ndt[d][t] + self.a) / (self.nd[d] + self.T * self.a)
    def estimate_phi(self):
        return [[self.pw_z(w, t) for w in range(self.W)] for t in range(self.T)]
    def estimate_theta(self):
        return [[self.pz_d(d, t) for t in range(self.T)] for d in range(self.D)]
    def topic_keys(self, num_displayed=5):
        phi, tks = self.estimate_phi(), []
        for w_ps in phi:
            tuples = [(p, self.vocab[i]) for i, p in enumerate(w_ps)]
            tuples.sort(reverse=True)
            tks.append([word for (p, word) in tuples[:num_displayed]])
        return tks
    def doc_keys(self, num_displayed=5, threshold=.02):
        theta, dks = self.estimate_theta(), []
        for t_ps in theta:
            tuples = [(p, t) for t, p in enumerate(t_ps)]
            tuples.sort(reverse=True)
            dks.append([(p, t) for (p, t) in tuples[:num_displayed] if p >= threshold])
        return dks
    def doc_detail(self, d):
        tks = self.topic_keys(num_displayed=5)
        for w, word in enumerate(self.docs[d]):
            topic, max_p = 0, 0
            for t in range(self.T):
                p = self.f(d, w, t)
                if p > max_p:
                    max_p, topic = p, t
            s = ' '.join(tks[topic])
            print '%s \t %i %s' % (word, topic, s)
def tokenize(docs):
    return [doc.split(" ") for doc in docs]

# docs = ['cat cat cat', 'cat dog rabbit']
docs = []
for line in fileinput.input("out.txt"):
    docs.append(line.strip())
fileinput.close()

tokenized_docs = tokenize(docs)
print tokenized_docs

lda = LDASampler(docs=tokenized_docs, num_topics=8, alpha=0.5, beta=0.5)
with open("wordmap.txt","w") as f:
    f.write(str(len(lda.wordmap()))+"\n")
    for k, v in lda.wordmap().iteritems():
        f.write(str(k)+" "+str(v)+"\n")

for i in xrange(500):
    print i
    lda.next()
with open("model-final.phi","w") as f:
    for item in lda.estimate_phi():
        f.write(" ".join([format(i,'.6f') for i in item])+"\n")
with open("model-final.theta","w") as f:
    for item in lda.estimate_theta():
        f.write(" ".join([format(i,'.6f') for i in item])+"\n")
