### DESCRIPTION

This is a toy example illustrating network streams in python.  We simulate data 
using the same toy procedure as in example1, then use a chime_packetizer transform 
to send it over the network, where it is received by a chime_network_stream object.
We put a waterfall_plotter in both the sending and receiving pipelines, so that
the plots can be compared as a visual test that the stream has been transmitted
correctly.

By default the sending and receiving pipelines are on the same node, and packets
are sent over the loopback address 127.0.0.1, but it's easy to modify the script
to run the pipelines on different machines.

### INSTRUCTIONS FOR RUNNING

```
    # Start the receiving pipeline first.  It will block, waiting for data...
    ./example5-receiver.py
    
    # ...in another window, start the sending pipeline.
    ./example5-sender.py
```

After both pipelines finish, you should see two sequences of waterfall plots 
sending_pipeline/waterfall_\*.png, receiving_pipeline/waterfall_\*.png 
which agree by eye.  Note that the simulated FRB appears in waterfall_5.png.

