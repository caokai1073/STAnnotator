"""Spatial cell type annotation engine with natural-language-driven constraints.

The LLM only translates user intent into tool calls; every annotation result is computed
deterministically here from (base probabilities, regions, constraints) and versioned.
Entry points: stagent.session.Session, stagent.mcp_server.
"""
__version__ = "0.2.42"  # one step up per release: 0.2.1, 0.2.2, ...
