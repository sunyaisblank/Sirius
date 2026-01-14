# Sirius Documentation

This directory contains the comprehensive documentation for the Sirius relativistic ray tracing engine. The documents are designed to be read in order, building from philosophical foundations through mathematical theory to practical operation.

## Document Overview

| Document | Purpose | Audience |
|----------|---------|----------|
| [Philosophy](philosophy.md) | Why the system exists; design principles | All readers |
| [Foundations](foundations.md) | Mathematical theory from first principles | Physicists, developers |
| [Types](types.md) | Type system, tensors, and domain modelling | Developers |
| [Specification](specification.md) | Formal system requirements | Implementers |
| [Guide](guide.md) | Practical operation | Operators |
| [Standard](standard.md) | Coding conventions | Developers |
| [Architecture](architecture.md) | System structure and naming | Developers |
| [Refactor](refactor.md) | Refactoring discipline | Developers |

## Reading Order

**For understanding the system:**

1. Philosophy: establishes the worldview and design rationale
2. Foundations: develops the mathematical machinery
3. Types: builds the type system from primitives to domain abstractions
4. Guide: explains practical operation

**For contributing to the codebase:**

1. Philosophy: understand the design principles
2. Types: understand the type system and invariants
3. Standard: learn the coding conventions
4. Architecture: understand the structural organisation

**For operating the renderer:**

1. Guide: primary reference for operation
2. Specification: reference for configuration and limits

## Document Summaries

### Philosophy

Explores the philosophical foundations: spacetime as geometry, the nature of light propagation, computational integrity, and the principles of high-integrity scientific visualisation. Reading this document illuminates why Sirius makes the design decisions it does.

### Foundations

Develops the mathematical theory from first principles: differential geometry, the geodesic equation, black hole spacetimes, exotic geometries, numerical integration, and radiative transfer. This document provides the theoretical grounding for all computational methods.

### Types

Defines the type system from primitives through domain abstractions, with guidance on coordinate representations, tensor types, and error handling patterns.

### Specification

Provides formal definitions of system behaviour: performance targets, mathematical invariants, precision requirements, test coverage mandates, and configuration parameters. This document serves as the authoritative reference for implementation correctness.

### Guide

Covers practical patterns for building, configuring, and operating the renderer. This document bridges theory and operation.

### Standard

Codifies conventions for high-integrity scientific software development: language compliance, predictable execution, defensive coding, GPU considerations, and testing requirements.

### Architecture

Describes the structural design and component responsibilities, including the component naming convention and registry.

### Refactor

Establishes the discipline for changing code whilst preserving integrity, performance, and auditability.
