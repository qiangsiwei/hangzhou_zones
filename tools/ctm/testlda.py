# -*- coding: utf-8 -*-

import itertools
import ldamodel, ctmmodel

def iterkeys(d):
	return iter(getattr(d, "iterkeys")())

def iteritems(d):
	return iter(getattr(d, "iteritems")())

class Dictionary(ldamodel.SaveLoad, dict):
	def __init__(self, documents=None):
		self.token2id = {}
		self.id2token = {}
		self.dfs = {}
		self.num_docs = 0
		self.num_pos = 0
		self.num_nnz = 0
		if documents is not None:
			self.add_documents(documents)
	def __getitem__(self, tokenid):
		if len(self.id2token) != len(self.token2id):
			self.id2token = dict((v, k) for k, v in iteritems(self.token2id))
		return self.id2token[tokenid]
	def keys(self):
		return self.token2id.values()
	def __len__(self):
		return len(self.token2id)
	def __str__(self):
		some_keys = list(itertools.islice(iterkeys(self.token2id), 5))
		return "Dictionary(%i unique tokens: %s%s)" % (len(self), some_keys, '...' if len(self) > 5 else '')
	def add_documents(self, documents):
		for docno, document in enumerate(documents):
			self.doc2bow(document, allow_update=True)
	def doc2bow(self, document, allow_update=False, return_missing=False):
		result = {}
		missing = {}
		document = sorted(token for token in document)
		for word_norm, group in itertools.groupby(document):
			frequency = len(list(group))
			tokenid = self.token2id.get(word_norm, None)
			if tokenid is None:
				if return_missing:
					missing[word_norm] = frequency
				if not allow_update:
					continue
				tokenid = len(self.token2id)
				self.token2id[word_norm] = tokenid
			result[tokenid] = frequency
		if allow_update:
			self.num_docs += 1
			self.num_pos += len(document)
			self.num_nnz += len(result)
			for tokenid in iterkeys(result):
				self.dfs[tokenid] = self.dfs.get(tokenid, 0) + 1
		result = sorted(iteritems(result))
		if return_missing:
			return result, missing
		else:
			return result

texts = [['human', 'interface', 'computer'], ['survey', 'user', 'computer', 'system', 'response', 'time'], ['eps', 'user', 'interface', 'system'], ['system', 'human', 'system', 'eps'], ['user', 'response', 'time'], ['trees'], ['graph', 'trees'], ['graph', 'minors', 'trees'], ['graph', 'minors', 'survey']]
dictionary = Dictionary(texts)
corpus = [dictionary.doc2bow(text) for text in texts]
# lda = ldamodel.LdaModel(corpus, num_topics=2, id2word=dictionary)
# print lda.print_topics()
lda = ctmmodel.CtmModel(corpus, num_topics=2, id2word=dictionary)
print lda.print_topics()

