tmhmmContainerFactory
========================
Usage
--------    
.. code-block:: python

   import pyproteinsExt.tmhmmContainer as tmhmm
   tmhmmContainer = tmhmm.parse('tmhmm.out')

You can iterate through tmhmmContainer

.. code-block:: python

    for e in tmhmmContainer:
        print(e.prot)
        print(e.prot_length)
        print(e.nb_helix)
        print(len(e.fragments),"fragments")
        for f in e.fragments:
            print(f.cellular_location, f.start, f.end)
        print(e.topology_seq)
        print()

Will print::

    tr|A0A2V3HRB7|A0A2V3HRB7_9EURY
    426
    0
    1 fragments
    outside 1 426
    iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

    tr|A0A2E0DFV0|A0A2E0DFV0_9EURY
    155
    4
    9 fragments
    outside 1 5
    TMhelix 6 25
    inside 26 53
    TMhelix 54 76
    outside 77 90
    TMhelix 91 109
    inside 110 129
    TMhelix 130 152
    outside 153 155
    iiiii11111111111111111111oooooooooooooooooooooooooooo22222222222222222222222iiiiiiiiiiiiii3333333333333333333oooooooooooooooooooo44444444444444444444444iii
    ...

You can access to one specific entry with its prot attribute

.. code-block:: python

   entry = tmhmmContainer.entries["tr|A0A2E0DFV0|A0A2E0DFV0_9EURY"]
   print(entry.prot, entry.prot_length, entry.nb_helix)

Will print::

    tr|A0A2E0DFV0|A0A2E0DFV0_9EURY 155 4

Content
--------
.. automodule:: tmhmmContainerFactory
    :members:
    :show-inheritance:
    :member-order: bysource