"""Spatial cell type annotation engine with natural-language-driven constraints.

The LLM only translates user intent into tool calls; every annotation result is computed
deterministically here from (base probabilities, regions, constraints) and versioned.
Entry points: celltypeagent.session.Session, celltypeagent.mcp_server.
"""
