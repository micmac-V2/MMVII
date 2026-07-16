::: center
**Help for command EditProfile**
:::

#  Descrition of EditProfile

This command inspects and edits MMVII user profiles. A user profile is a
set of key/value preferences used by MMVII commands, for example the
default serialization formats, the maximum number of parallel processes,
or the command used to open PDF files.

When EditProfile is called without option, it prints:

-   the known profile names, with the current profile shown between
    square brackets;

-   the path of the profile file effectively loaded;

-   the key/value pairs of the current profile.

```{=html}
<!-- -->
```
    MMVII EditProfile

# Selecting the current profile

The optional parameter `SetCurrent` makes a profile the default profile
for subsequent MMVII commands:

    MMVII EditProfile SetCurrent=Default
    MMVII EditProfile SetCurrent=Test

The profile name must contain only letters, digits, underscores, or
hyphens. It must not be empty. If the selected profile file does not
exist yet, MMVII creates it by synchronizing it with the distributed
default profile.

The current profile name is stored in the user configuration file
`MMVII-Current-Profile.xml`.

# Editing profile values

The optional parameter `KeyVal` sets one key/value pair in the selected
profile. It must contain exactly two values: the key and the value. To
specify an empty string for a value, use the special syntax
`KeyVal=[Key,,]`.

    MMVII EditProfile KeyVal=[NbProcMax,8]
    MMVII EditProfile KeyVal=[TaggedSerialMode,json]
    MMVII EditProfile KeyVal=[VectSerialMode,csv]
    MMVII EditProfile KeyVal=[PdfOpen,evince]
    MMVII EditProfile KeyVal=[Empty,,]

The command writes the modified profile file after applying the change.
Values are stored as strings in XML and converted later by the code that
reads the key.

The optional parameter `DelKey` removes a key from the selected profile:

    MMVII EditProfile DelKey=PdfOpen

Removing a key from a user profile does not remove it from the
distributed default profile. At the next synchronization, keys present
in `MMVII-Default-Profile.xml` are restored with their default values
unless the user profile overrides them.

# Command options

::: center
  ----------------------------------------------------------------
  Option         Type         Meaning
  -------------- ------------ ------------------------------------
  `SetCurrent`   profile name Name of the profile to make current

  `KeyVal`       vector of    Set a profile value as `[Key,Value]`
                 two strings  

  `DelKey`       profile key  Delete one key from the selected
                              profile
  ----------------------------------------------------------------
:::

If none of these options is used, EditProfile only displays the current
profile information and does not modify any profile file.

# Using a profile for one command

Every MMVII command has a global `Profile` option. It selects a profile
only for the current execution:

    MMVII SomeCommand ... Profile=Test

This does not modify `MMVII-Current-Profile.xml`. It is useful for
scripts or for testing profile-specific behavior without changing the
user's default profile.

Command-line options remain authoritative for one command execution. For
example, `NbProcMax` is only applied when the command did not explicitly
initialize the global `NbProc` option. When both are present, the
explicit command-line option keeps priority.

# User profile mechanism

MMVII user profiles store per-user defaults outside the project data
directory. They are loaded by every MMVII command. A command normally
uses the current profile, but it can also select another profile for one
execution with the global `Profile` option.

## Configuration directory

User profiles are stored in an operating-system-specific MMVII
configuration directory:

::: center
  -------------------------------------------------------------
  System    Directory
  --------- ---------------------------------------------------
  Linux     `$XDG_CONFIG_HOME/MMVII` when `XDG_CONFIG_HOME` is
            an absolute path, otherwise `$HOME/.config/MMVII`

  macOS     `$HOME/Library/Application Support/MMVII`

  Windows   `%APPDATA%/MMVII`
  -------------------------------------------------------------
:::

MMVII creates this directory when it needs to write the current profile
selection or a user profile file.

## Profile files

The profile mechanism uses three file families:

  ----------------------------------------------------------------
  File                          Location
  ----------------------------- ----------------------------------
  `MMVII-Default-Profile.xml`   `MMVII-LocalParameters/`
                                installation directory

  `MMVII-Current-Profile.xml`   user configuration directory

  `MMVII-profile-<Name>.xml`    user configuration directory
  ----------------------------------------------------------------

  : Configuration files and their locations

  ----------------------------------------------------------------
  File                          Role
  ----------------------------- ----------------------------------
  `MMVII-Default-Profile.xml`   Default key/value set distributed
                                with MMVII

  `MMVII-Current-Profile.xml`   Stores the name of the current
                                profile.

  `MMVII-profile-<Name>.xml`    Stores the user-defined values for
                                profile `<Name>`.
  ----------------------------------------------------------------

  : Configuration files and their roles

If `MMVII-Current-Profile.xml` does not exist, MMVII uses the profile
name `Default`. When a profile is loaded, MMVII first reads
`MMVII-Default-Profile.xml`, then overlays values from
`MMVII-profile-<Name>.xml` if that file already exists. It writes the
merged result back to the user profile file. This keeps user profile
files synchronized when new default keys are added by MMVII.

## Standard keys

The default profile currently defines these keys:

::: center
  ---------------------------------------------------------------------
  Key                  Default   Meaning
  -------------------- --------- --------------------------------------
  `NbProcMax`          `4`       Upper limit for the number of parallel
                                 processes used by commands when
                                 `NbProc` was not explicitly provided

  `TaggedSerialMode`   `xml`     Default serialization mode for tagged
                                 objects

  `VectSerialMode`     `csv`     Default serialization mode for
                                 vector-like serialized data

  `PdfOpen`            empty     Command used by MMVII to open
                       string    generated PDF files, when configured
  ---------------------------------------------------------------------
:::

If a key is absent, MMVII uses the caller-provided default value and
emits a user warning. The profile is initialized early enough for help
and argument specification generation, so profile names and profile keys
can be exposed as completion or GUI metadata through argument types
tagged as `eTA2007::Profile` and `eTA2007::ProfileKey`.
