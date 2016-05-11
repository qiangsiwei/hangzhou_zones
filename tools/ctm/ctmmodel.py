# -*- coding: utf-8 -*-

import copy
import numpy
from numpy.linalg import inv, det
from scipy.optimize import minimize, fmin_l_bfgs_b

class SufficientStats():
    def __init__(self, numtopics, numterms):
        self.numdocs = 0
        self.numtopics = numtopics
        self.numterms = numterms
        self.beta_stats = numpy.zeros([numtopics, numterms])
        self.mu_stats = numpy.zeros(numtopics)
        self.sigma_stats = numpy.zeros([numtopics, numtopics])
    def update(self, lamda, nu2, phi, doc):
        """
        Given optimized variational parameters, update statistics
        """
        # update mu_stats
        self.mu_stats += lamda
        # update \beta_stats[i], 0 < i < self.numtopics
        for n, c in doc:
            for i in xrange(self.numtopics):
                self.beta_stats[i, n] += c * phi[n, i]
        # update \sigma_stats
        self.sigma_stats += numpy.diag(nu2) + numpy.dot(lamda, lamda.transpose())
        self.numdocs += 1

class CtmModel():
    def __init__(self, corpus=None, num_topics=100, id2word=None, estep_convergence=0.001, em_convergence=0.0001, em_max_iterations=10):
        self.id2word = id2word
        self.estep_convergence = estep_convergence
        self.em_convergence = em_convergence
        self.em_max_iterations = em_max_iterations
        self.num_terms = 1 + max(self.id2word.keys())
        self.num_topics = int(num_topics)
        # initialize a model with zero-mean, diagonal covariance gaussian and random topics seeded from the corpus
        self.mu = numpy.zeros(self.num_topics)
        self.sigma = numpy.diagflat([1.0] * self.num_topics)
        self.sigma_inverse = inv(self.sigma)
        self.beta = numpy.random.uniform(0, 1, (self.num_topics, self.num_terms))
        # variational parameters
        self.lamda = numpy.zeros(self.num_topics)
        self.nu2 = numpy.ones(self.num_topics)  # nu^2
        self.phi = 1/float(self.num_topics) * numpy.ones([self.num_terms, self.num_topics])
        self.optimize_zeta()
        # in order to get the topics graph, we need to store the optimized lamda for each document
        self.observed_lamda = numpy.zeros([len(corpus)])
        # if a training corpus was provided, start estimating the model right away
        if corpus is not None:
            self.expectation_maximization(corpus)

    def expectation_maximization(self, corpus):
        """
        Expectation-Maximization algorithm.
        During E-step, variational parameters are optimized with fixed model parameters.
        During M-step, model parameters are optimized given statistics collected in E-step.
        """
        for iteration in xrange(self.em_max_iterations):
            old_bound = self.corpus_bound(corpus)
            print iteration
            # print "bound before E-step %f" %(old_bound)
            # E-step and collect sufficient statistics for the M-step
            statistics = self.do_estep(corpus)
            # M-step
            self.do_mstep(statistics)
            new_bound = self.corpus_bound(corpus)
            # print "bound after M-step %f" %(new_bound)
            if (new_bound - old_bound)/old_bound < self.em_convergence:
                break

    def do_estep(self, corpus):
        # initialize empty statistics
        statistics = SufficientStats(self.num_topics, self.num_terms)
        for d, doc in enumerate(corpus):
            # variational_inference modifies the variational parameters
            model = copy.deepcopy(self)
            model.variational_inference(doc)
            # collect statistics for M-step
            statistics.update(model.lamda, model.nu2, model.phi, doc)
        return statistics

    def do_mstep(self, sstats):
        """
        Optimize model's parameters using the statictics collected during the e-step
        """
        for i in xrange(self.num_topics):
            beta_norm = numpy.sum(sstats.beta_stats[i])
            self.beta[i] = sstats.beta_stats[i] / beta_norm
        self.mu = sstats.mu_stats / sstats.numdocs
        self.sigma = sstats.sigma_stats + numpy.multiply(self.mu, self.mu.transpose())
        self.sigma_inverse = inv(self.sigma)

    def bound(self, doc, lamda=None, nu2=None):
        """
        Estimate the variational bound of a document
        """
        if lamda is None:lamda = self.lamda
        if nu2 is None:nu2 = self.nu2
        N = sum([cnt for _, cnt in doc])  # nb of words in document
        bound = 0.0
        # E[log p(\eta | \mu, \Sigma)] + H(q(\eta | \lamda, \nu) + sum_n,i { \phi_{n,i}*log(\phi_{n,i}) }
        bound += 0.5 * numpy.log(det(self.sigma_inverse))
        bound -= 0.5 * numpy.trace(numpy.dot(numpy.diag(nu2), self.sigma_inverse))
        bound -= 0.5 * (lamda - self.mu).transpose().dot(self.sigma_inverse).dot(lamda - self.mu)
        bound += 0.5 * (numpy.sum(numpy.log(nu2)) + self.num_topics) # TODO safe_log
        # print "first term %f for doc %s" %(bound, doc)
        # \sum_n { E[log p(z_n | \eta)] - sum_i {\lamda_i * \phi_{n, i}}
        sum_exp = numpy.sum([numpy.exp(lamda[i] + 0.5*nu2[i]) for i in xrange(self.num_topics)])
        bound += (N * (-1/self.zeta * sum_exp + 1 - numpy.log(self.zeta)))
        # print "second term %f for doc %s" %(bound, doc)
        # E[log p(w_n | z_n, \beta)] - sum_n,i { \phi_{n,i}*log(\phi_{n,i}) }
        bound += sum([c * self.phi[n, i] * (lamda[i] + numpy.log(self.beta[i, n]) - numpy.log(self.phi[n, i]))
            for (n, c) in doc
            for i in xrange(self.num_topics)
        ])
        return bound

    def corpus_bound(self, corpus):
        """
        Estimates the likelihood bound for the whole corpus by summing over all the documents in the corpus.
        """
        return sum([self.bound(doc) for doc in corpus])

    def variational_inference(self, doc):
        """
        Optimize variational parameters (zeta, lamda, nu, phi) given the current model and a document. This method modifies the model self.
        """
        bound = self.bound(doc)
        new_bound = bound
        # print "bound before variational inference %f" %(bound)
        for iteration in xrange(self.em_max_iterations):
            # print "bound before zeta opt %f" %(self.bound(doc))
            self.optimize_zeta()
            # print "bound after zeta opt %f" %(self.bound(doc))
            # print "bound before lamda opt %f" %(self.bound(doc))
            self.optimize_lamda(doc)
            # print "bound after lamda opt %f" %(self.bound(doc))
            self.optimize_zeta()
            # print "bound before nu2 opt %f" %(self.bound(doc))
            self.optimize_nu2(doc)
            # print "bound after nu2 opt %f" %(self.bound(doc))
            self.optimize_zeta()
            # print "bound before phi opt %f" %(self.bound(doc))
            self.optimize_phi(doc)
            # print "bound after phi opt %f" %(self.bound(doc))
            bound, new_bound = new_bound, self.bound(doc)
            relative_change = abs((bound - new_bound)/bound)
            if (relative_change < self.estep_convergence):
                break
        # print "bound after variational inference %f" %(bound)
        return bound

    def optimize_zeta(self):
        self.zeta = sum([numpy.exp(self.lamda[i] + 0.5 * self.nu2[i]) for i in xrange(self.num_topics)])

    def optimize_phi(self, doc):
        for n, _ in doc:
            phi_norm = sum([numpy.exp(self.lamda[i]) * self.beta[i, n] for i in xrange(self.num_topics)])
            for i in xrange(self.num_topics):
                self.phi[n, i] = numpy.exp(self.lamda[i]) * self.beta[i, n] / phi_norm

    def optimize_lamda(self, doc):
        def f(lamda):
            return self.bound(doc, lamda=lamda)
        def df(lamda):
            """
            Returns dL/dlamda
            """
            N = sum([c for _, c in doc])
            result = numpy.zeros(self.num_topics)
            result -= numpy.dot(self.sigma_inverse, (lamda - self.mu))
            result += sum([c * self.phi[n, :] for n, c in doc])
            result -= (N/self.zeta)*numpy.array([numpy.exp(lamda[i] + 0.5 * self.nu2[i]) for i in xrange(self.num_topics)])
            return result
        # We want to maximize f, but numpy only implements minimize, so we minimize -f
        res = minimize(lambda x: -f(x), self.lamda, method='BFGS', jac=lambda x: -df(x))
        self.lamda = res.x

    def optimize_nu2(self, doc):
        def f(nu2):
            return self.bound(doc, nu2=nu2)
        def df(nu2):
            """
            Returns dL/dnu2
            """
            N = sum([c for _, c in doc])
            result = numpy.zeros(self.num_topics)
            for i in xrange(self.num_topics):
                result[i] = - 0.5 * self.sigma_inverse[i, i]
                result[i] -= N/(2*self.zeta) * numpy.exp(self.lamda[i] + 0.5 * nu2[i])
                result[i] += 1/(2*nu2[i]) # TODO safe_division
            return result
        # constraints : we need nu2[i] >= 0
        constraints = [(0, None) for _ in xrange(self.num_topics)]
        result = fmin_l_bfgs_b(lambda x: -f(x), self.nu2, fprime=lambda x: -df(x), bounds=constraints)
        self.nu2 = result[0]

    def print_topics(self):
        print "self.mu", self.mu
        print "self.sigma", self.sigma
        print "self.beta", self.beta
        print "self.lamda", self.lamda
        print "self.nu2", self.nu2
        print "self.phi", self.phi
        print " self.observed_lamda", self.observed_lamda
