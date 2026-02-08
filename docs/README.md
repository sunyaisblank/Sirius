# Documentation

Documentation for the Sirius relativistic ray tracing engine. The documents are ordered so that each builds on the previous.

## Documents

| Document | What it covers |
|----------|----------------|
| [Design Rationale](philosophy.md) | Why the system is built the way it is: physical assumptions, computational strategy, design principles |
| [Foundations](foundations.md) | Mathematical theory: differential geometry, geodesic equation, spacetimes, numerical methods |
| [Types](type.md) | Type system: numeric primitives, coordinates, tensors, metric interface, ray state |
| [Specification](specification.md) | Performance targets, invariant tolerances, precision requirements, test coverage |
| [Guide](guide.md) | Building, configuring, and running the renderer |
| [Standard](standard.md) | Coding conventions (JPL/MISRA-based, C++17) |
| [Architecture](architecture.md) | Layer structure, naming convention, component registry |
| [Refactoring](refactor.md) | Method for changing code while preserving correctness and performance |
| [Commentary](comment.md) | Inline documentation and comment standards |

## Reading Paths

**Understanding the system:** Design Rationale, then Foundations, then Types, then Guide.

**Contributing code:** Design Rationale, then Types, then Standard, then Architecture.

**Operating the renderer:** Guide, then Specification for configuration limits.
