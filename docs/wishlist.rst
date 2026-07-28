Wishlist
========

isEE_framework is a cleaned up and simplified version of the root isEE package. It's still not publishable as-is. Here, I am documenting a short list of changes that will need to be made before a version of isEE_framework can be published as part of a standalone package.

Probably essential changes

* Remove pytraj dependency. Versions of pytraj distributed through conda and pypi are extremely restrictive on compatible Python versions. Direct installation is a big ask and a major pain point. Most usages of pytraj in isEE_framework ought to be replaced with calls to MDTraj instead.
* Replace hard-coded tleap template with optional user-provided template.
* Resolve all or almost all outstanding todos, especially those denoting kludges
* The strategy for implementing mutations (which is obviously highly central to isEE_framework in any form) is deeply ugly and fragile. Consider removing everything related to ts_bonds
* Fix the entrypoint script so that isEE_framework can be invoked directly as a command, rather than with ``python path/to/main.py``

Valuable but maybe not essential changes

* Remove internal dependency on Amber as the MD engine. Everything that currently depends on Amber-specific file formats or functions should be moved to methods of the Amber implementation of the MDEngine abstract base class in infrastructre/mdengine. At a minimum we should also support Gromacs.
* The entire initialize_charges.py file in particular is highly Amber-specific (and deeply sloppy anyway). It could be removed; I have no idea if there's a legitimate use-case fot it.
* Remove dependency on pAPRika; find another way to hook into tleap (if tleap is still being used internally). Maybe call it on the command line with subprocess. pAPRika is a nightmare for the dependency chain and seriously limits portability.
* stabilitymodel.py and related settings probably go in the trash.
* Figure out whether there's a cleaner way to depend on pyrosetta, since it's license is controlled