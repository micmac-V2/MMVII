import sys
import os
import json
import subprocess
import re
import io
import contextlib
import builtins
import sys

cur_word=sys.argv[1]
comp_cword=int(sys.argv[2])
#  Optional: a directory holding the output of the two MMVII calls this script would
#  otherwise make, under the names read by readCacheFile below. It is what lets the
#  completion work where MMVII cannot be run, typically a host completing an MMVII
#  that only exists inside a container.
cache_dir=sys.argv[3] if len(sys.argv)>3 else None
comp_line=os.getenv('COMP_LINE')
screen_columns=int(os.getenv('COLUMNS','80'))

# For windows, force remove CR (^M) from output
sys.stdout.reconfigure(newline="\n")

helpArgPossible=True
help_printed=False
completion_prefix=""
completion_nospace=False
completion_suffix=""
structured_completion=False

debug=0
debug_tty=os.getenv('MMVII_COMPLETION_DEBUG_TTY')
if debug_tty:
    sys.stderr = open(debug_tty, "w", buffering=1)
    debug=1


if debug:
    _original_print = builtins.print

    def print(*args, **kwargs):
        # Affichage normal sur stdout
        _original_print(*args, **kwargs)
        # Même affichage sur stderr
        kwargs2 = kwargs.copy()
        kwargs2["file"] = sys.stderr
        _original_print(*args, **kwargs2)

    def printd(*args, **kwargs):
        kwargs2 = kwargs.copy()
        kwargs2["file"] = sys.stderr
        _original_print(*args, **kwargs2)

    def dbg(*args, **kwargs):
        kwargs2 = kwargs.copy()
        kwargs2["file"] = sys.stderr
        _original_print("DEBUG <<<",file=sys.stderr)
        _original_print(*args, **kwargs2)
        _original_print(">>>DEBUG",file=sys.stderr)

def printCompletionOptions(options='') -> None:
    if structured_completion:
        options = options.replace('-o filenames','')
    if completion_nospace:
        options += ' -o nospace'
    options += ' -o nosort'
    #  compopt takes '-o name' pairs and repeating one is meaningless, while the same option
    #  legitimately comes from two places : an ambiguous value adds '-o nospace' that its
    #  caller may already have given. Keep each name once, in the order it first appeared
    names = re.findall(r'-o\s+(\S+)',options)
    print ("Options:" + ' '.join('-o ' + name for name in dict.fromkeys(names)))


def completedValue(value) -> str:
    return completion_prefix + value + completion_suffix


def printCompletedValue(word, matches, hide_brackets=False, suffixes=None) -> bool:
    """Print the completion candidates, return True when the value stays unfinished.

       Readline replaces the whole word by the longest common prefix of the candidates, so a
       candidate that does not begin with what is already typed truncates the line. They all
       carry completion_prefix for that reason. When the common prefix adds something, one
       single candidate is enough : it is inserted, and nothing is listed. A help line already
       printed forbids that shortcut : being a candidate of its own, it empties the prefix
       readline computes, so nothing would be inserted and the values would just be hidden.
    """
    if len(matches) == 1:
        match = re.sub(r'^=','',matches[0])
        #  suffixes : the separator that follows depends on the value itself, as for the name of
        #  an interpolator, which is followed by ',' or by the closing of the expression
        if suffixes and (match in suffixes):
            print(completion_prefix + match + suffixes[match])
        else:
            print(completedValue(match))
        return False

    common = os.path.commonprefix(matches)
    if not help_printed and len(common) > len(word):
        print(completion_prefix + common)
        return True

    for w in matches:
        if hide_brackets:
            w = w.rstrip('[')
        print(re.sub(r'^=','',completion_prefix + word + w[len(word):]))
    return True


def hasUnclosedQuote(value) -> bool:
    quote = None
    escaped = False
    for char in value or '':
        if escaped:
            escaped = False
        elif char == '\\':
            escaped = True
        elif quote:
            if char == quote:
                quote = None
        elif char in "'\"":
            quote = char
    return quote is not None


def printFinishedValue(value,quoted_only=False) -> bool:
    if not completion_suffix or not value:
        return False
    if quoted_only and not (comp_line or "").rstrip().endswith(("'","\"")):
        return False
    print(completedValue(value))
    printCompletionOptions()
    return True


#  The interpolators declared by MMVII, read from the "config" section of GenArgsSpec.
#  An argument having the "Interpol" semantic is a flat prefixed expression of this language,
#  so its completion is driven by these specifications and never by a duplicated list
theInterpolatorList = []
theInterpolators = {}


def isInterpolSpec(spec) -> bool:
    return 'Interpol' in (spec.get('semantic') or [])


def listElements(value) -> list:
    """The elements of a bracketed list, the last one being the one under edition"""
    if not value.startswith('['):
        return []
    result = []
    depth = 0
    begin = 1
    end = len(value)
    for index,char in enumerate(value[1:],1):
        if char in '[{(':
            depth += 1
        elif char in ']})':
            if depth == 0:
                end = index
                break
            depth -= 1
        elif char == ',' and depth == 0:
            result.append(value[begin:index])
            begin = index + 1
    result.append(value[begin:end])
    return result


def interpolRoles(elements) -> tuple:
    """Role of each position of the expression, given the elements already typed : ('name',) or
       ('param',interpolator,parameter). Same walk as cInterpolator1D::AllocFromNames : a name,
       its parameters, then the interpolator it wraps. Stops as soon as a name is missing or
       unknown; the second result tells whether the expression is closed at that length."""
    roles = []
    pending = 1
    while pending > 0:
        position = len(roles)
        roles.append(('name',))
        if position >= len(elements):
            return roles,False
        interpol = theInterpolators.get(elements[position])
        if interpol is None:
            return roles,False
        pending -= 1
        for param in interpol['params']:
            roles.append(('param',interpol,param))
        if interpol['subInterpol']:
            pending += 1
    return roles,True


def interpolElementSpec(role) -> dict:
    if role[0] == 'name':
        #  the names themselves are the completion values, no need to describe them one by one
        return {'name':'Interpolator','type':'string','comment':'name of the interpolator',
                'allowed':[interpol['name'] for interpol in theInterpolatorList],
                'interpolElement':True}
    _,interpol,param = role
    return {'name':interpol['name'] + '.' + param['name'],'type':param['type'],
            'interpolElement':True,
            'comment':param['comment'] + ' (e.g. ' + param['example'] + ')'}


def isVectorSpec(spec) -> bool:
    return bool(spec.get('isVector')) or re.fullmatch(r'(?:std::)?vector<.*>',spec.get('type','')) is not None


def vectorElementSpec(spec) -> dict:
    result = dict(spec)
    result['isVector'] = False
    result.pop('vsize',None)
    vector = re.fullmatch(r'(?:std::)?vector<(.*)>',spec.get('type',''))
    if vector:
        result['type'] = vector.group(1)
        # The help line must keep showing the type of the field, which is a vector : only the
        # value being completed is one of its elements
        result['fieldType'] = '[' + vector.group(1) + ',...]'
    return result


def currentListElement(value) -> tuple:
    if not value.startswith('['):
        return 'invalid',
    depth = 0
    element = 0
    begin = 1
    for index,char in enumerate(value[1:],1):
        if char in '[{(':
            depth += 1
        elif char in ']})':
            if depth == 0:
                if index == len(value)-1:
                    return 'complete',
                return 'invalid',
            depth -= 1
        elif char == ',' and depth == 0:
            element += 1
            begin = index + 1
    return 'active',element,value[begin:],value[:begin]


def structuredBracketDepth(spec) -> int:
    depth = 0
    element_spec = spec
    while True:
        if isVectorSpec(element_spec):
            depth += 1
            element_spec = vectorElementSpec(element_spec)
            continue
        fields = element_spec.get('fields')
        if fields:
            depth += 1
            element_spec = fields[0]
            continue
        break
    return depth


def structuredOpeningPrefix(spec) -> str:
    return '[' * max(0,structuredBracketDepth(spec)-1)


def structuredFieldOpening(spec) -> str:
    # Brackets opening a structured field, empty for a plain one. They are inserted together
    # with the ',' that reaches the field, so that it does not take one completion more
    return '[' * structuredBracketDepth(spec)


def structuredFieldContext(spec,value,prefix='',top_level=True,close_suffix='') -> dict:
    # A top-level vector argument with CanRepeat can be given as a singleton,
    # i.e. with one less leading '[' than the fully-bracketed vector form (see
    # bracket_depth/V_InitParam in MMVII_Tpl_ElemStrToVal.h). When that applies,
    # delegate entirely to the element type, one bracket level down.
    if top_level and isVectorSpec(spec) and '##CanRepeat' in (spec.get('semantic') or []):
        leading = len(value) - len(value.lstrip('['))
        if not value or leading == structuredBracketDepth(spec)-1:
            element_spec = vectorElementSpec(spec)
            element_spec.pop('semantic',None)
            return structuredFieldContext(element_spec,value,prefix,True)
    fields = spec.get('fields')
    if not fields and not isVectorSpec(spec):
        return {'kind':'field','spec':spec,'value':value,'prefix':prefix,'suffix':'','finish':top_level}
    if not value:
        return {'kind':'opening','prefix':prefix+structuredOpeningPrefix(spec)}

    current = currentListElement(value)
    if current[0] == 'complete':
        return {'kind':'complete','value':prefix+value}
    if current[0] != 'active':
        return {'kind':'invalid'}
    _,index,field_value,value_prefix = current
    prefix += value_prefix

    if isVectorSpec(spec):
        max_elements = parseList(spec.get('vsize'))[1]
        can_repeat = max_elements is None or index+1 < max_elements
        element_spec = vectorElementSpec(spec)
        completing_name = False
        must_repeat = False
        if isInterpolSpec(spec):
            elements = listElements(value)
            roles,_ = interpolRoles(elements[:index])
            if index >= len(roles):        # nothing more is expected by the grammar
                return {'kind':'invalid'}
            element_spec = interpolElementSpec(roles[index])
            completing_name = roles[index][0] == 'name'
            if completing_name:
                #  what follows a name depends on the name itself : ',' when the grammar expects
                #  more, else the closing of the expression and of the level that contains it
                element_spec['allowedSuffix'] = {
                    name : (',' if len(interpolRoles(elements[:index]+[name])[0]) > index+1
                                else ']' + close_suffix)
                    for name in theInterpolators
                }
            #  once the current element is known, the grammar says whether another one is
            #  required (only ','), or whether the expression is closed (only ']')
            roles,closed = interpolRoles(elements[:index] + [field_value])
            must_repeat = len(roles) > index+1
            can_repeat = must_repeat or not (closed and len(roles) == index+1)
        child = structuredFieldContext(element_spec,field_value,prefix,False)
        if child['kind'] == 'complete':
            return {'kind':'hold','value':child['value'],'can_repeat':can_repeat,
                    'must_repeat':must_repeat,'finish':top_level,'close_suffix':close_suffix}
        if completing_name and (field_value not in theInterpolators):
            return child                   # let the name be completed before offering a separator
        if not element_spec.get('fields') and not isVectorSpec(element_spec) and child['kind'] == 'field' and child['value']:
            return {'kind':'hold','value':child['prefix']+child['value'],'can_repeat':can_repeat,
                    'must_repeat':must_repeat,'finish':top_level,'close_suffix':close_suffix}
        return child
    if index >= len(fields):
        return {'kind':'invalid'}

    suffix = ',' if index+1 < len(fields) else ']'
    finish = top_level and suffix == ']'
    if suffix == ',':
        suffix += structuredFieldOpening(fields[index+1])
    field = fields[index]
    # A field is structured when it is a nested struct or a vector : both have their own
    # bracket level, so they must be handled by a recursive call and not as a plain value
    if field.get('fields') or isVectorSpec(field):
        child = structuredFieldContext(field,field_value,prefix,False,suffix)
        if child['kind'] == 'complete':
            return {'kind':'append','value':child['value'],'suffix':suffix,'finish':finish}
        if (child['kind'] == 'hold') and (not child['can_repeat']):
            child['finish'] = finish       # closing the field also closes or separates the level
        group = {'fields':fields,'owner':spec.get('name')}
        child['field_groups'] = [group] + child.get('field_groups',[])
        return child
    return {
        'kind':'field',
        'spec':field,
        'value':field_value,
        'prefix':prefix,
        'suffix':suffix,
        'finish':finish,
        'field_groups':[{'fields':fields,'owner':spec.get('name')}]
    }

def parseList(s) -> list:
    if not s:
        return [None,None]
    return [int(x) if x!='' else None for x in s[1:-1].split(',')]

def printLineWide(s) -> None:
    if not s:
        return
    global help_printed
    help_printed = True
    if isinstance(s, str):
        s = [s]
    for line in s:
        spaces_needed = screen_columns - len(line) % screen_columns - 2
        print (line + ' ' * spaces_needed)

def print_line(s) -> None:
    if not s:
        return
    printLineWide(s)
    print(' ')

def print_help_line(s) -> None:
    printLineWide(s)

def fieldNameType(field) -> str:
    name = field.get('name','')
    optional = name.endswith('?')
    if optional:
        name = name[:-1]
    return f"{name}:{field.get('type','')}" + ('?' if optional else '')


def structuredHeader(spec,field_groups) -> list:
    comment = spec.get('comment','structured argument')
    name = spec.get('name')
    if isVectorSpec(spec):
        msg = [f"> Vector of structure : {comment}"]
    elif name:
        msg = [f"> Structure for [{name}]: {comment}"]
    else:
        msg = [f"> Structure expected : {comment}"]
    for index,group in enumerate(field_groups):
        names = ', '.join(fieldNameType(field) for field in group['fields'])
        owner = '' if index == 0 else group.get('owner')
        suffix = f" (for {owner})" if owner else ''
        msg.append(f"> Fields: [{names}]{suffix}")
    return msg

def printMsgExit(s,code=0) -> None:
    print_line(s)
    printCompletionOptions()
    sys.exit(code)


def printHelps(word) -> bool:
    if not helpArgPossible or len(word)==0:
        return False
    matches=[w for w in ('help','Help','HELP') if w.startswith(word)]
    if len(matches) == 1:
        print(completedValue(re.sub(r'^=','',matches[0])))
        return True
    return False
    

def printFilter(word,words,option='',hide_brackets=False,suffixes=None) -> None:
    word_lower = word.lower()
    matches=[w for w in sorted(words) if w.lower().startswith(word_lower)]
    if len(matches) == 0:
        printHelps(word)
        printCompletionOptions()
        return
    else:
        if printCompletedValue(word,matches,hide_brackets,suffixes):
            option += " -o nospace"    # an ambiguous value is never finished
    if suffixes and len(matches)==1:
        #  the value is finished when its own suffix closes it, whatever the spec said
        global completion_nospace
        completion_nospace = not suffixes.get(matches[0],'').endswith(']')
    printCompletionOptions(option)


def printFiles(word,extensions=None,only_dirs=False) -> None:
    word = os.path.expanduser(word)
    (path,name)=os.path.split(word)
    search_path = path if path != '' else '.'

    dirs  = []
    files = []
    try:
        with os.scandir(search_path) as it:
            for file in it:
                if file.name.startswith('.') and not name.startswith('.'):
                    continue
                try:
                    if file.is_dir():
                        dirs.append(os.path.join(path,file.name))
                    elif not only_dirs and (not extensions or "" in extensions or any (file.name.lower().endswith(ext) for ext in extensions)):
                        files.append(os.path.join(path,file.name))
                except PermissionError:
                    pass
    except (FileNotFoundError,PermissionError):
        pass
    matches=[w for w in sorted(files)+sorted(dirs) if w.startswith(word)]
    option = '-o filenames'
    if len(matches) == 0:
        if printFinishedValue(word,True):
            return
        printHelps(word)
        printCompletionOptions()
        return
    else:
        if printCompletedValue(word,matches):
            option += " -o nospace"    # an ambiguous value is never finished
    printCompletionOptions(option)


def printDirs(word,path) -> None:
    word = os.path.expanduser(word)
    dirs = []
    try:
        with os.scandir(path) as it:
            for file in it:
                if file.name.startswith('.') and not word.startswith('.'):
                    continue
                try:
                    if file.is_dir():
                        dirs.append(file.name)
                except PermissionError:
                    pass
    except (FileNotFoundError,PermissionError):
        pass

    matches=[w for w in sorted(dirs) if w.startswith(word)]
    option = '-o filenames'
    if len(matches) == 0:
        if printFinishedValue(word,True):
            return
        printHelps(word)
        printCompletionOptions()
        return
    else:
        if printCompletedValue(word,matches):
            option += " -o nospace"    # an ambiguous value is never finished
    printCompletionOptions(option)



def printSpecHelp(all_specs, spec, value) -> None:
    global completion_prefix,completion_nospace,completion_suffix,structured_completion
    msg=[]
    vector_header = False
    is_struct_field = False
    if spec.get('fields') or isVectorSpec(spec):
        root_spec = spec
        structured_completion = True
        context = structuredFieldContext(spec,value)
        kind = context['kind']
        if kind == 'opening':
            completion_prefix = context['prefix']
            completion_nospace = True
            printFilter('',['['],'-o nospace')
            return
        if kind == 'field':
            spec = context['spec']
            value = context['value']
            completion_prefix = context['prefix']
            completion_suffix = context['suffix']
            completion_nospace = not context['finish']
            is_struct_field = bool(root_spec.get('fields'))
            if not value and root_spec.get('fields'):
                msg += structuredHeader(root_spec,context.get('field_groups',[]))
            elif not value and isVectorSpec(root_spec) and not spec.get('interpolElement'):
                msg.append(f'>Expect <{root_spec["type"].replace("std::","")}>: {root_spec.get("comment","vector argument")}')
                vector_header = True
        elif kind == 'append':
            print(context['value']+context['suffix'])
            options = '' if context['finish'] else '-o nospace'
            printCompletionOptions(options)
            return
        elif kind == 'hold':
            if context.get('must_repeat'):
                #  only ',' is possible : print it alone, a help line would be one more completion
                #  candidate and bash would list them instead of inserting the separator
                print(context['value']+',')
                printCompletionOptions('-o nospace')
                return
            if context['can_repeat']:
                print_help_line(">Expect <separator>: ',' for another value or ']' to close the vector")
                print(context['value']+',')
            else:
                print(context['value']+']'+context.get('close_suffix',''))
                printCompletionOptions('' if context['finish'] else '-o nospace')
                return
            print(context['value']+']')
            printCompletionOptions('-o nospace')
            return
        elif kind == 'complete':
            print(context['value'])
            printCompletionOptions()
            return
        else:
            return
    atype = spec['type']
    msg_type = atype.replace('std::','')
    semantic = list(spec.get('semantic') or [])
    allowed = spec.get('allowed')
    vrange = parseList(spec.get('range'))
    vsize = parseList(spec.get('vsize'))
    vector = re.search(r'vector<(.*)>',msg_type)
    if vector:
        msg_type='[' + vector.group(1) + ',...]'
    elif atype == 'bool':
        allowed=['True','1','False','0']

    if not vector_header:
        if is_struct_field:
            field_name = spec.get('name','')
            optional = field_name.endswith('?')
            if optional:
                field_name = field_name[:-1]
            msg.append(f'> Current field: {field_name}: {spec.get("fieldType",msg_type)}' + ('?' if optional else '') +
                       (f' : {spec.get("comment")}' if spec.get("comment") else '') +
                       (f' [Default: {spec.get("default")}]' if spec.get("default") else ''))
        else:
            # A mandatory (positional) argument has no 'name': fall back to its comment
            # instead of a filler string, as a struct field's own name never comes up here.
            name = spec.get("name")
            comment = spec.get("comment")
            desc = f'{name} : {comment}' if name and comment else name or comment or ''
            msg.append(f'>Expect <{msg_type}>: {desc}' +
                       (f' [Default: {spec.get("default")}]' if spec.get("default") else ''))

    if vector and vsize[0] and vsize[0]>1:
        if not value:
            printFilter (value,'[','-o nospace')
            return
        printMsgExit(msg)

    if not value:
        print_line(msg)

    if vrange[0] and vrange[1]:
        if vrange[1] - vrange[0] < 20:
            printFilter(value,[str(i) for i in range(vrange[0],vrange[1]+1)])
        return

    if allowed:
        printFilter(value,allowed,suffixes=spec.get('allowedSuffix'))
        return
        
    if msg_type != 'string':
        if printFinishedValue(value):
            return
        if not printHelps(value):
            printMsgExit('')
        return
    if not semantic:
        if printFinishedValue(value,True):
            return
        if not printHelps(value):
            printMsgExit('')
        return

    dir_types = all_specs['config']['eTa2007DirTypes']
    dir_type = set(semantic).intersection(dir_types)
    if len(dir_type) == 1:
        dir_type = list(dir_type)[0]
        sub_dir = all_specs['config']['MMVIIDirPhp'] +  dir_type
        printDirs(value,sub_dir)
        return

    file_types = all_specs['config']['eTa2007FileTypes']
    if 'MPF' in semantic:
        semantic.append('Im')
    file_types = set(semantic).intersection(file_types)
    if len(file_types):
        extensions = []
        for file_type in file_types:
            extensions += all_specs['config']['extensions'][file_type]
        if 'MPF' in semantic:
            extensions += ['.xml','.json']
        printFiles(value,extensions)
        return
    file_types = ('FDP','In','Out','OptEx')
    if len(set(semantic).intersection(file_types)):
        printFiles(value)
        return
    if "DP" in semantic:
        if completion_prefix:
            printFiles(value,only_dirs=True)
        else:
            print(f"File:compgen -d -o filenames -- '{value}'")
        return

def readCacheFile(name) -> str:
    """ Content of a file of cache_dir, standing for the output of an MMVII call """
    try:
        with open(os.path.join(cache_dir,name)) as f:
            return f.read()
    except OSError:
        printMsgExit(f">ERROR: can't read {name} in {cache_dir}.",1)


def printPatBench(arg) -> None:
    if cache_dir:
        output=readCacheFile('BenchNames.txt')
    else:
        try:
            result=subprocess.run(['MMVII','Bench','1','PatBench=XXX'],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
        except:
            printMsgExit('>ERROR: MMVII not found.',1)
        output=result.stdout

    benches=re.findall(r'^\s*[-#]\s*(\S+)', output, flags=re.MULTILINE)
    benches.append("XXX")
    printFilter(arg,benches)


def getAllSpecs() -> dict:
    #  Holds the command list too, under 'applets': there is no separate source for it
    if cache_dir:
        output=readCacheFile('GenArgsSpec.json')
    else:
        try:
            result=subprocess.run(['MMVII','GenArgsSpec'],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
        except:
            printMsgExit('>ERROR: MMVII not found.',1)
        output=result.stdout
    try:
        all_specs=json.loads(output)
    except json.JSONDecodeError as e:
        if debug:
            print("JSON error: ",e,sys.stderr)
        printMsgExit(">ERROR: Can't get args specification from MMVII.",1)
    return all_specs

def commandNames(applets) -> None:
    app_names=(a['name'] for a in applets)
    printFilter(cur_word,app_names)

def optionalArgCompletion(spec) -> str:
    result = spec['name'] + '='
    if spec.get('fields') or isVectorSpec(spec):
        context = structuredFieldContext(spec,'')
        result += context['prefix']
        if context['kind'] == 'opening':
            result += '['
    return result

def main() -> int:
    if hasUnclosedQuote(comp_line):
        return 0
    all_specs=getAllSpecs()
    global theInterpolatorList,theInterpolators
    theInterpolatorList = (all_specs.get('config') or {}).get('interpolators') or []
    theInterpolators = {interpol['name']:interpol for interpol in theInterpolatorList}
    try:
        applets=all_specs['applets']
    except:
        printMsgExit(">ERROR: Can't get args specification from MMVII.",1)
    
    # 1st argument: MMVII command name
    if comp_cword == 1:
        commandNames(applets)
        return 0

    # else get applet specification
    command=comp_line.split()[1].lower()
    applet=[a for a in applets if a['name'].lower() == command]
    if len(applet) != 1:
        sys.exit(0)
    applet=applet[0]

    # Nth first arguments are mandatory
    if comp_cword < len(applet['mandatory'])+2:
        spec = applet['mandatory'][comp_cword-2]
        printSpecHelp(all_specs, spec, cur_word)
        return 0

    # Others are optionals
    arg_split=re.search(r'([-+a-zA-Z0-9_.]+)(=(.*))?',cur_word)
    #  no '='
    if arg_split==None or arg_split.start(2) < 0:
        normal=sorted((optionalArgCompletion(a) for a in applet['optional'] if a['level'] == 'normal'),key=str.lower)
        tuning=sorted((optionalArgCompletion(a) for a in applet['optional'] if a['level'] == 'tuning'),key=str.lower)
        common=sorted((optionalArgCompletion(a) for a in applet['optional'] if a['level'] == 'global'),key=str.lower)
        printFilter(cur_word,normal+tuning+common,'-o nospace',hide_brackets=True)
        return 0
    # got '='
    global helpArgPossible
    helpArgPossible = False
    arg=arg_split.group(1)
    if arg == 'PatBench' or arg == 'PatRefutBench=':
        printPatBench(arg_split.group(3))
        return 0
    specs = [ s for s in applet.get('optional') if s['name'].lower() == arg.lower() ] 
    if len(specs) != 1:
        return 0
    printSpecHelp(all_specs, specs[0], arg_split.group(3))
    return 0



if __name__ == '__main__':
    sys.exit(main())

