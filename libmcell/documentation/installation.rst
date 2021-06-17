*******************
MCell4 Installation
*******************

Download
########

CellBlender bundle containig MCell4 is available for download 
for different operating systems at the MCell.org website: `<https://mcell.org/download.html>`_.

Installation
############

MacOS
*****

After downloading, the zip file will be automatically extracted to your downloads directory. 
Move the Blender-2.79-CellBlender directory to Applications (i.e. to the directory /Applications). 
The CellBlender bundle won’t work correctly in any other directory.

MCell4 build is known to work also on the latest Apple MacBooks with ARM processor M1.

MacOS Mojave and Older
----------------------
 
If you have an older version of MacOS than the Catalina MacOS, you can skip  
to the following section `Running Blender`_ because the following setup is not needed and the only 
thing required is that the unpacked bundle is under the directory /Applications.
 

Start Blender by double-clicking the blender application file. 

.. image:: images/installation_macos_start_blender.png

Although the package is digitally signed, some newer MacOS versions require extra approvals from Apple, 
for now, you will most probably get one of the following warnings:

.. image:: images/installation_macos_warning1.png

.. image:: images/installation_macos_warning2.png


Click on **Cancel** or on **OK** and then open System Preferences (through the Apple menu in the top left). 
Select Security and Privacy.

.. image:: images/installation_macos_security_and_privacy.png

In the Security & Privacy settings click on **Open Anyway**.

.. image:: images/installation_macos_security_and_privacy_selected.png

One more warning appears, select **Open**.

.. image:: images/installation_macos_warning3.png

.. image:: images/installation_macos_warning4.png

Now quit Blender and start it from the terminal as shown in the following section. 
This will allow you to see additional messages printed by CellBlender. 
It also allows one to overcome a settings saving issue on MacOS Catalina (described in section 2.2).


Running Blender
---------------

Open a terminal window, terminal can be found under Applications and Utilities.

.. image:: images/installation_macos_start_terminal.png


Run the following commands from the terminal:

.. code-block:: text

      cd /
      /Applications/Blender-2.79-CellBlender/my_blender

By now, CellBlender should be up and running, however, if you get a message that the application 
is damaged, please see section `Common Troubleshooting`_.


Setting System Variable MCELL_PATH
##################################

TODO


Common Troubleshooting
######################

TODO


