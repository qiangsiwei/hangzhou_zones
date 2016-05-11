# -*- coding: utf-8 -*-

import itertools

import numpy
from scipy.special import gammaln, digamma, psi
from scipy.special import gamma as gammafunc
from scipy.special import polygamma
try:
    from scipy.maxentropy import logsumexp
except ImportError:
    from scipy.misc import logsumexp

def grouper(iterable, chunksize, as_numpy=False):
    it = iter(iterable)
    while True:
        if as_numpy:
            wrapped_chunk = [[numpy.array(doc) for doc in itertools.islice(it, int(chunksize))]]
        else:
            wrapped_chunk = [list(itertools.islice(it, int(chunksize)))]
        if not wrapped_chunk[0]:
            break
        yield wrapped_chunk.pop()

def is_corpus(obj):
    try:
        if 'Corpus' in obj.__class__.__name__:
            return True, obj
    except:
        pass
    try:
        if hasattr(obj, 'next'):
            doc1 = next(obj)
            obj = itertools.chain([doc1], obj)
        else:
            doc1 = next(iter(obj))
        if len(doc1) == 0:
            return True, obj
        id1, val1 = next(iter(doc1))
        id1, val1 = int(id1), float(val1)
    except:
        return False, obj
    return True, obj

def dirichlet_expectation(alpha):
    if (len(alpha.shape) == 1):
        result = psi(alpha) - psi(numpy.sum(alpha))
    else:
        result = psi(alpha) - psi(numpy.sum(alpha, 1))[:, numpy.newaxis]
    return result.astype(alpha.dtype)

class SaveLoad(object):
    @classmethod
    def load(cls, fname, mmap=None):
        subname = lambda suffix: fname + '.' + suffix + '.npy'
        obj = unpickle(fname)
        for attrib in getattr(obj, '__numpys', []):
            setattr(obj, attrib, numpy.load(subname(attrib), mmap_mode=mmap))
        for attrib in getattr(obj, '__scipys', []):
            sparse = unpickle(subname(attrib))
            sparse.data = numpy.load(subname(attrib) + '.data.npy', mmap_mode=mmap)
            sparse.indptr = numpy.load(subname(attrib) + '.indptr.npy', mmap_mode=mmap)
            sparse.indices = numpy.load(subname(attrib) + '.indices.npy', mmap_mode=mmap)
            setattr(obj, attrib, sparse)
        for attrib in getattr(obj, '__ignoreds', []):
            setattr(obj, attrib, None)
        return obj
    def save(self, fname, separately=None, sep_limit=10 * 1024**2, ignore=frozenset()):
        subname = lambda suffix: fname + '.' + suffix + '.npy'
        tmp = {}
        if separately is None:
            separately = []
            for attrib, val in iteritems(self.__dict__):
                if isinstance(val, numpy.ndarray) and val.size >= sep_limit:
                    separately.append(attrib)
                elif isinstance(val, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)) and val.nnz >= sep_limit:
                    separately.append(attrib)
        for attrib in separately + list(ignore):
            if hasattr(self, attrib):
                tmp[attrib] = getattr(self, attrib)
                delattr(self, attrib)
        try:
            numpys, scipys, ignoreds = [], [], []
            for attrib, val in iteritems(tmp):
                if isinstance(val, numpy.ndarray) and attrib not in ignore:
                    numpys.append(attrib)
                    numpy.save(subname(attrib), numpy.ascontiguousarray(val))
                elif isinstance(val, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)) and attrib not in ignore:
                    scipys.append(attrib)
                    numpy.save(subname(attrib) + '.data.npy', val.data)
                    numpy.save(subname(attrib) + '.indptr.npy', val.indptr)
                    numpy.save(subname(attrib) + '.indices.npy', val.indices)
                    data, indptr, indices = val.data, val.indptr, val.indices
                    val.data, val.indptr, val.indices = None, None, None
                    try:
                        pickle(val, subname(attrib))
                    finally:
                        val.data, val.indptr, val.indices = data, indptr, indices
                else:
                    ignoreds.append(attrib)
            self.__dict__['__numpys'] = numpys
            self.__dict__['__scipys'] = scipys
            self.__dict__['__ignoreds'] = ignoreds
            pickle(self, fname)
        finally:
            for attrib, val in iteritems(tmp):
                setattr(self, attrib, val)

class LdaState(SaveLoad):
    def __init__(self, eta, shape):
        self.eta = eta
        self.sstats = numpy.zeros(shape)
        self.numdocs = 0
    def reset(self):
        self.sstats[:] = 0.0
        self.numdocs = 0
    def merge(self, other):
        assert other is not None
        self.sstats += other.sstats
        self.numdocs += other.numdocs
    def blend(self, rhot, other, targetsize=None):
        assert other is not None
        if targetsize is None:
            targetsize = self.numdocs
        if self.numdocs == 0 or targetsize == self.numdocs:
            scale = 1.0
        else:
            scale = 1.0 * targetsize / self.numdocs
        self.sstats *= (1.0 - rhot) * scale
        if other.numdocs == 0 or targetsize == other.numdocs:
            scale = 1.0
        else:
            scale = 1.0 * targetsize / other.numdocs
        self.sstats += rhot * scale * other.sstats
        self.numdocs = targetsize
    def blend2(self, rhot, other, targetsize=None):
        assert other is not None
        if targetsize is None:
            targetsize = self.numdocs
        self.sstats += other.sstats
        self.numdocs = targetsize
    def get_lambda(self):
        return self.eta + self.sstats
    def get_Elogbeta(self):
        return dirichlet_expectation(self.get_lambda())

class CorpusABC(SaveLoad):
    def save(self, *args, **kwargs):
        super(CorpusABC, self).save(*args, **kwargs)
    @staticmethod
    def save_corpus(fname, corpus, id2word=None, metadata=False):
        with open(fname, 'w') as fout:
            for doc in corpus:
                fmt = str(doc)
                fout.write("%s\n" % fmt)

class TransformedCorpus(CorpusABC):
    def __init__(self, obj, corpus, chunksize=None):
        self.obj, self.corpus, self.chunksize = obj, corpus, chunksize
        self.metadata = False
    def __len__(self):
        return len(self.corpus)
    def __iter__(self):
        if self.chunksize:
            for chunk in grouper(self.corpus, self.chunksize):
                for transformed in self.obj.__getitem__(chunk, chunksize=None):
                    yield transformed
        else:
            for doc in self.corpus:
                yield self.obj[doc]

class TransformationABC(SaveLoad):
    def _apply(self, corpus, chunksize=None):
        return TransformedCorpus(self, corpus, chunksize)

class LdaModel(TransformationABC):
    def __init__(self, corpus=None, num_topics=100, id2word=None, distributed=False, chunksize=2000, passes=1, update_every=1, alpha='symmetric', eta=None, decay=0.5, eval_every=10, iterations=50, gamma_threshold=0.001):
        self.id2word = id2word
        self.num_terms = 1 + max(self.id2word.keys())
        self.distributed = bool(distributed)
        self.num_topics = int(num_topics)
        self.chunksize = chunksize
        self.decay = decay
        self.num_updates = 0
        self.passes = passes
        self.update_every = update_every
        self.eval_every = eval_every
        self.optimize_alpha = alpha == 'auto'
        self.alpha = numpy.asarray([1.0 / num_topics for i in xrange(num_topics)])
        if eta is None:
            self.eta = 1.0 / num_topics
        else:
            self.eta = eta
        self.iterations = iterations
        self.gamma_threshold = gamma_threshold
        self.dispatcher = None
        self.numworkers = 1
        self.state = LdaState(self.eta, (self.num_topics, self.num_terms))
        self.state.sstats = numpy.random.gamma(100., 1. / 100., (self.num_topics, self.num_terms))
        self.sync_state()
        if corpus is not None:
            self.update(corpus)

    def sync_state(self):
        self.expElogbeta = numpy.exp(self.state.get_Elogbeta())

    def clear(self):
        self.state = None
        self.Elogbeta = None

    def inference(self, chunk, collect_sstats=False):
        try:
            _ = len(chunk)
        except:
            chunk = list(chunk)
        gamma = numpy.random.gamma(100., 1. / 100., (len(chunk), self.num_topics))
        Elogtheta = dirichlet_expectation(gamma)
        expElogtheta = numpy.exp(Elogtheta)
        if collect_sstats:
            sstats = numpy.zeros_like(self.expElogbeta)
        else:
            sstats = None
        converged = 0
        for d, doc in enumerate(chunk):
            ids = [id for id, _ in doc]
            cts = numpy.array([cnt for _, cnt in doc])
            gammad = gamma[d, :]
            Elogthetad = Elogtheta[d, :]
            expElogthetad = expElogtheta[d, :]
            expElogbetad = self.expElogbeta[:, ids]
            phinorm = numpy.dot(expElogthetad, expElogbetad) + 1e-100
            for _ in xrange(self.iterations):
                lastgamma = gammad
                gammad = self.alpha + expElogthetad * numpy.dot(cts / phinorm, expElogbetad.T)
                Elogthetad = dirichlet_expectation(gammad)
                expElogthetad = numpy.exp(Elogthetad)
                phinorm = numpy.dot(expElogthetad, expElogbetad) + 1e-100
                meanchange = numpy.mean(abs(gammad - lastgamma))
                if (meanchange < self.gamma_threshold):
                    converged += 1
                    break
            gamma[d, :] = gammad
            if collect_sstats:
                sstats[:, ids] += numpy.outer(expElogthetad.T, cts / phinorm)
        if collect_sstats:
            sstats *= self.expElogbeta
        return gamma, sstats

    def do_estep(self, chunk, state=None):
        if state is None:
            state = self.state
        gamma, sstats = self.inference(chunk, collect_sstats=True)
        state.sstats += sstats
        state.numdocs += gamma.shape[0]
        return gamma

    def update_alpha(self, gammat, rho):
        N = float(len(gammat))
        logphat = sum(dirichlet_expectation(gamma) for gamma in gammat) / N
        dalpha = numpy.copy(self.alpha)
        gradf = N * (psi(numpy.sum(self.alpha)) - psi(self.alpha) + logphat)
        c = N * polygamma(1, numpy.sum(self.alpha))
        q = -N * polygamma(1, self.alpha)
        b = numpy.sum(gradf / q) / ( 1 / c + numpy.sum(1 / q))
        dalpha = -(gradf - b) / q
        if all(rho() * dalpha + self.alpha > 0):
            self.alpha += rho() * dalpha
        return self.alpha

    def log_perplexity(self, chunk, total_docs=None):
        if total_docs is None:
            total_docs = len(chunk)
        corpus_words = sum(cnt for document in chunk for _, cnt in document)
        subsample_ratio = 1.0 * total_docs / len(chunk)
        perwordbound = self.bound(chunk, subsample_ratio=subsample_ratio) / (subsample_ratio * corpus_words)
        return perwordbound

    def update(self, corpus, chunksize=None, decay=None, passes=None, update_every=None, eval_every=None, iterations=None, gamma_threshold=None):
        chunksize = self.chunksize
        decay = self.decay
        passes = self.passes
        update_every = self.update_every
        eval_every = self.eval_every
        iterations = self.iterations
        gamma_threshold = self.gamma_threshold
        rho = lambda: pow(1.0 + self.num_updates, -decay)
        lencorpus = len(corpus)
        self.state.numdocs += lencorpus
        updatetype = "batch"
        updateafter = lencorpus
        evalafter = min(lencorpus, (eval_every or 0) * self.numworkers * chunksize)
        updates_per_pass = max(1, lencorpus / updateafter)
        for pass_ in xrange(passes):
            if self.dispatcher:
                self.dispatcher.reset(self.state)
            else:
                other = LdaState(self.eta, self.state.sstats.shape)
            dirty = False
            reallen = 0
            for chunk_no, chunk in enumerate(grouper(corpus, chunksize, as_numpy=True)):
                reallen += len(chunk)
                if eval_every and ((reallen == lencorpus) or ((chunk_no + 1) % (eval_every * self.numworkers) == 0)):
                    self.log_perplexity(chunk, total_docs=lencorpus)
                if self.dispatcher:
                    self.dispatcher.putjob(chunk)
                else:
                    gammat = self.do_estep(chunk, other)
                    if self.optimize_alpha:
                        self.update_alpha(gammat, rho)
                dirty = True
                del chunk
                if update_every and (chunk_no + 1) % (update_every * self.numworkers) == 0:
                    if self.dispatcher:
                        other = self.dispatcher.getstate()
                    self.do_mstep(rho(), other)
                    del other
                    if self.dispatcher:
                        self.dispatcher.reset(self.state)
                    else:
                        other = LdaState(self.eta, self.state.sstats.shape)
                    dirty = False
            if dirty:
                if self.dispatcher:
                    other = self.dispatcher.getstate()
                self.do_mstep(rho(), other)
                del other
                dirty = False

    def do_mstep(self, rho, other):
        diff = numpy.log(self.expElogbeta)
        self.state.blend(rho, other)
        del other
        diff -= self.state.get_Elogbeta()
        self.sync_state()
        self.print_topics(15)
        self.num_updates += 1

    def bound(self, corpus, gamma=None, subsample_ratio=1.0):
        score = 0.0
        _lambda = self.state.get_lambda()
        Elogbeta = dirichlet_expectation(_lambda)
        for d, doc in enumerate(corpus):
            if gamma is None:
                gammad, _ = self.inference([doc])
            else:
                gammad = gamma[d]
            Elogthetad = dirichlet_expectation(gammad)
            score += numpy.sum(cnt * logsumexp(Elogthetad + Elogbeta[:, id]) for id, cnt in doc)
            score += numpy.sum((self.alpha - gammad) * Elogthetad)
            score += numpy.sum(gammaln(gammad) - gammaln(self.alpha))
            score += gammaln(numpy.sum(self.alpha)) - gammaln(numpy.sum(gammad))
        score *= subsample_ratio
        score += numpy.sum((self.eta - _lambda) * Elogbeta)
        score += numpy.sum(gammaln(_lambda) - gammaln(self.eta))
        score += numpy.sum(gammaln(self.eta * self.num_terms) - gammaln(numpy.sum(_lambda, 1)))
        return score

    def print_topics(self, topics=10, topn=10):
        return self.show_topics(topics, topn, log=True)

    def show_topics(self, topics=10, topn=10, log=False, formatted=True):
        if topics < 0 or topics >= self.num_topics:
            topics = self.num_topics
            chosen_topics = range(topics)
        else:
            topics = min(topics, self.num_topics)
            sort_alpha = self.alpha + 0.0001 * numpy.random.rand(len(self.alpha))
            sorted_topics = list(numpy.argsort(sort_alpha))
            chosen_topics = sorted_topics[ : topics/2] + sorted_topics[-topics/2 : ]
        shown = []
        for i in chosen_topics:
            if formatted:
                topic = self.print_topic(i, topn=topn)
            else:
                topic = self.show_topic(i, topn=topn)
            shown.append(topic)
        return shown

    def show_topic(self, topicid, topn=10): 
        topic = self.state.get_lambda()[topicid]
        topic = topic / topic.sum()
        bestn = numpy.argsort(topic)[::-1][:topn]
        beststr = [(topic[id], self.id2word[id]) for id in bestn]
        return beststr

    def print_topic(self, topicid, topn=10):
        return ' + '.join(['%.3f*%s' % v for v in self.show_topic(topicid, topn)])

    def __getitem__(self, bow, eps=0.01):
        is_corpus, corpus = is_corpus(bow)
        if is_corpus:
            return self._apply(corpus)
        gamma, _ = self.inference([bow])
        topic_dist = gamma[0] / sum(gamma[0])
        return [(topicid, topicvalue) for topicid, topicvalue in enumerate(topic_dist) if topicvalue >= eps]

    def save(self, fname, *args, **kwargs):
        if self.state is not None:
            self.state.save(fname + '.state', *args, **kwargs)
        super(LdaModel, self).save(fname, *args, ignore=['state', 'dispatcher'], **kwargs)

    @classmethod
    def load(cls, fname, *args, **kwargs):
        kwargs['mmap'] = kwargs.get('mmap', 'r')
        result = super(LdaModel, cls).load(fname, *args, **kwargs)
        try:
            result.state = super(LdaModel, cls).load(fname + '.state', *args, **kwargs)
        except Exception as e:
            logging.warning("failed to load state from %s: %s" % (fname + '.state', e))
        return result
