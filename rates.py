from abc import ABCMeta, abstractmethod
from collections import namedtuple
import math


RateParams = namedtuple('RateParams', ['rate', 'scale', 'midpoint'])


class RateFunctions:

	def __init__(self, mathlib):
		self.math = math

	def explin(self, v, params, tol=1e-3):
		a, k, vo = params
		δv = v - vo
		if abs(δv) < tol:
			z = a * k
		else:
			v1 = k * δv
			if v1 > 0:
				exp_v = self.math.exp(-v1)
				z = δv * exp_v / (1 - exp_v)
			else:
				z = δv / (self.math.exp(v1) - 1)
		return z * a

	def _sig_n(self, x):
		return 1 / (1 + self.math.exp(x))

	def _sig_p(self, x):
		exp_x = self.math.exp(-x)
		return exp_x / (1 + exp_x)

	def logistic(self, v, params):
		m = self.math
		a, k, vo = params
		δv = (v - vo) / k
		z = m.where(m.abs(δv) < 0,
		            self._sig_n(δv),
		            self._sig_p(δv))
		return z * a

	def exp(self, v, params):
		a, k, vo = params
		δv = (v - vo) / k
		return a * self.math.exp(δv)


class Rate(metaclass=ABCMeta):

	def __init__(self, params: RateParams):
		self.params = params

	@abstractmethod
	def __call__(self, v):
		pass


